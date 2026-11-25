"""Run one case and coupling resolution in its own AMUSE process."""
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import json
import os
from pathlib import Path
import re
import traceback

import numpy as np

from survey_orbits import (
    FH_WD_TYPES, STELLAR_STRUCTURE_FIELDS, fh_coupling_grid, fh_crop_track,
    fh_driver, fh_json, fh_run_case,
)


@contextmanager
def case_lock(path):
    """Hold an OS-released advisory lock, including across an interrupted run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as lock:
        if os.name == "nt":
            import msvcrt
            if path.stat().st_size == 0:
                lock.write(b"0")
                lock.flush()
            lock.seek(0)
            try:
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as exc:
                raise RuntimeError(f"A writer already holds {path}") from exc
        else:
            import fcntl
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError(f"A writer already holds {path}") from exc
        yield


def load_driver(path, evolving_tide=False, pilot_age_Myr=None):
    fields = ["time_yr", "mass_MSun", "radius_AU", "stellar_type", "age_Myr"]
    if evolving_tide:
        fields.extend(STELLAR_STRUCTURE_FIELDS)
    with np.load(path, allow_pickle=False) as source:
        missing = set(fields) - set(source.files)
        if missing:
            raise ValueError(f"Stellar driver lacks {sorted(missing)}. Run prepare for this tide model.")
        track = fh_driver({field: source[field] for field in fields})
    if track["age_Myr"][0] != 0:
        raise ValueError("The survival experiment starts at stellar age zero.")
    if pilot_age_Myr is not None:
        if not 0 < pilot_age_Myr <= track["age_Myr"][-1]:
            raise ValueError("The pilot endpoint must be positive and within the stellar driver.")
        if pilot_age_Myr < track["age_Myr"][-1]:
            track = fh_crop_track(track, end_age_Myr=pilot_age_Myr)
    else:
        wd = np.isin(track["stellar_type"], FH_WD_TYPES)
        first_wd = int(np.argmax(wd))
        if first_wd == 0 or first_wd == len(wd) - 1 or not wd[first_wd:].all():
            raise ValueError("Full runs need a progenitor followed by a continuous WD phase "
                             "and at least one post-transition sample.")
    return track


def _prepare_case(folder, metadata):
    record = folder / "job.json"
    if record.is_file():
        if json.loads(record.read_text(encoding="utf-8")) != metadata:
            raise ValueError(f"Settings differ from {folder}. Choose a new output for changed settings.")
        if not (folder / "checkpoint.json").is_file():
            status_file = folder / "run.json"
            if status_file.is_file():
                status = json.loads(status_file.read_text(encoding="utf-8")).get("status")
                if status == "completed":
                    raise RuntimeError(f"Missing terminal checkpoint for {folder}")
            # A crash during worker startup can occur before the initial checkpoint.
            attempts = folder.parent / ".attempts"
            attempts.mkdir(exist_ok=True)
            number = 1
            while (attempts / f"{folder.name}-{number}").exists():
                number += 1
            folder.rename(attempts / f"{folder.name}-{number}")
    elif folder.is_dir() and any(folder.iterdir()):
        raise RuntimeError(f"Unrecognized files in existing case directory {folder}")

    folder.mkdir(parents=True, exist_ok=True)
    if not record.is_file():
        fh_json(record, metadata)


def run_job(case_dir, driver_path, spec, settings, level, history_stride, *,
            max_outputs=25_000_000, pilot_age_Myr=None, suite=None):
    """Return the checked summary. Identical invocations reuse or resume a case."""
    match = re.fullmatch(r"n([1-9][0-9]*)(_tight)?", level)
    if not match:
        raise ValueError(f"Invalid coupling resolution: {level!r}")
    subdivisions = int(match.group(1))
    factor = settings.get("tight_factor", 10.0) if match.group(2) else 1.0
    folder = Path(case_dir).resolve()
    metadata = dict(schema="planet_survival_job_v1", case=spec["label"],
                    data_level=level, spec=spec, suite=suite or spec.get("suite"),
                    history_stride=history_stride, pilot_age_Myr=pilot_age_Myr, settings=settings)
    with case_lock(folder.parent / ".locks" / (folder.name + ".lock")):
        driver = load_driver(driver_path, spec["cfg"].get("stellar_tide_prescription") == 1,
                             pilot_age_Myr=pilot_age_Myr)
        targets = fh_coupling_grid(
            driver, spec["step_yr"], subdivisions=subdivisions, max_outputs=max_outputs,
            max_mass_fraction=settings.get("max_mass_fraction", 0.002),
            max_radius_fraction=settings.get("max_radius_fraction", 0.01),
        )
        _prepare_case(folder, metadata)
        with (folder / "worker.log").open("a", buffering=1, encoding="utf-8") as log:
            with redirect_stdout(log), redirect_stderr(log):
                try:
                    _, summary = fh_run_case(
                        folder, driver, spec["cfg"], targets,
                        {**settings["physics"], "companion": spec["inclination_deg"] is not None},
                        rtol=settings["rtol"] / factor, atol=settings["atol"] / factor,
                        worker_options=settings["worker_options"], limits=settings["limits"],
                        resume=True, history_stride=history_stride,
                        pilot=pilot_age_Myr is not None,
                    )
                except BaseException:
                    traceback.print_exc()
                    raise
    return summary
