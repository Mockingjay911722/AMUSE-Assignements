"""Build the shared SeBa stellar driver and the envelope inputs for Hurley tides."""

from array import array
from contextlib import chdir
from importlib.metadata import PackageNotFoundError, version
import json
import math
import os
from pathlib import Path
import platform
import tempfile
from time import perf_counter

import numpy as np

from orbit_worker import case_lock
from survey_orbits import FH_RSUN_AU, FH_WD_TYPES, STELLAR_STRUCTURE_FIELDS, fh_driver

DRIVER_FIELDS = ("time_yr", "mass_MSun", "radius_AU", "stellar_type", "age_Myr")
BASE_STEPS = {"0": 5., "1": 5., "2": 0.2, "3": 0.2, "4": 0.5, "5": 0.01, "6": 0.01}
DEFAULTS = dict(initial_mass_MSun=1.2, metallicity=0.02, pre_agb_refinement=128,
                agb_refinement=2048, wd_step_yr=0.5, wd_after_yr=100.,
                max_age_Myr=14_000., base_steps_Myr=BASE_STEPS)
SCHEMA = "mixed_seba_driver_v1"


def normalized_settings(overrides):
    """Return validated settings in a form that also round-trips through JSON."""
    unknown = set(overrides) - DEFAULTS.keys()
    if unknown:
        raise ValueError(f"Unknown stellar settings: {sorted(unknown)}")
    settings = {**DEFAULTS, **overrides}
    for key in ("pre_agb_refinement", "agb_refinement"):
        value = settings[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{key} must be a positive integer")
    steps = {str(k): v for k, v in settings["base_steps_Myr"].items()}
    if set(steps) != set(BASE_STEPS):
        raise ValueError("base_steps_Myr needs entries for SeBa types 0 through 6")
    for key, value in list(settings.items()) + list(steps.items()):
        if key != "base_steps_Myr" and (
            isinstance(value, bool) or not isinstance(value, (float, int))
            or not math.isfinite(value) or value <= 0
        ):
            raise ValueError(f"{key} must be a finite positive number")
    return {**settings, "base_steps_Myr": steps}


def _post_wd_requests(first_wd, duration, step):
    # Anchor to the first WD age. Repeated floating-point addition drifts across the tail.
    end = first_wd + duration / 1e6
    if end <= first_wd:
        raise ValueError("wd_after_yr is too small to resolve at this stellar age")
    for i in range(1, math.ceil(duration / step)):
        target = first_wd + i * step / 1e6
        if target < end:
            yield target
    yield end


def _generate(settings, channel, until_age_Myr):
    from amuse.community.seba import Seba
    from amuse.datamodel import Particle
    from amuse.units import units

    # Columns of doubles keep memory bounded even for the 1.67-million-row reference.
    fields = ("age_Myr", "mass_MSun", "radius_AU", "stellar_type",
              *STELLAR_STRUCTURE_FIELDS, "requested_time_Myr")
    columns = {field: array("d") for field in fields}
    steps = {int(k): v for k, v in settings["base_steps_Myr"].items()}
    last_report = started = perf_counter()
    stellar = Seba(channel_type=channel)
    try:
        stellar.parameters.metallicity = settings["metallicity"]
        star = stellar.particles.add_particle(Particle(mass=settings["initial_mass_MSun"] | units.MSun))

        def record(target):
            nonlocal last_report
            age = float(stellar.model_time.value_in(units.Myr))
            kind = float(star.stellar_type.value_in(units.stellar_type))
            if (not math.isfinite(age) or abs(age - target) > 1e-9 or not kind.is_integer()
                    or (columns["age_Myr"] and age <= columns["age_Myr"][-1])):
                raise RuntimeError(f"Unexpected SeBa time or phase at requested age {target} Myr")
            values = (age, star.mass.value_in(units.MSun), star.radius.value_in(units.RSun) * FH_RSUN_AU,
                      kind, star.luminosity.value_in(units.LSun),
                      star.convective_envelope_mass.value_in(units.MSun),
                      star.convective_envelope_radius.value_in(units.RSun), target)
            for field, value in zip(fields, values):
                columns[field].append(value)
            if perf_counter() - last_report >= 30:
                print(f"SeBa: {len(columns['age_Myr']):,} states, stellar age {age:.6f} Myr, "
                      f"{(perf_counter() - started) / 60:.1f} min elapsed", flush=True)
                last_report = perf_counter()
            return age, int(kind)

        age, kind = record(0.)
        refined = False
        while kind not in FH_WD_TYPES and age != until_age_Myr:
            if kind not in steps or age >= settings["max_age_Myr"]:
                raise RuntimeError(f"SeBa reached stellar type {kind} at {age} Myr before WD formation")
            refined = refined or kind in (5, 6)
            level = settings["agb_refinement"] if refined else settings["pre_agb_refinement"]
            target = min(age + steps[kind] / level, settings["max_age_Myr"])
            if until_age_Myr is not None:
                target = min(target, until_age_Myr)
            stellar.evolve_model(target | units.Myr)
            age, kind = record(target)
            if target == until_age_Myr:
                break
        if kind in FH_WD_TYPES and (until_age_Myr is None or age < until_age_Myr):
            duration = settings["wd_after_yr"]
            if until_age_Myr is not None and until_age_Myr > age + duration / 1e6:
                raise ValueError("The pilot age exceeds the configured white-dwarf endpoint")
            wd_targets = _post_wd_requests(age, duration, settings["wd_step_yr"])
        else:
            wd_targets = ()
        for target in wd_targets:
            if until_age_Myr is not None:
                target = min(target, until_age_Myr)
            stellar.evolve_model(target | units.Myr)
            _, kind = record(target)
            if kind not in FH_WD_TYPES:
                raise RuntimeError(f"SeBa left the WD phase at {target} Myr")
            if target == until_age_Myr:
                break
    finally:
        stellar.stop()

    result = {k: np.frombuffer(v, dtype=np.float64) for k, v in columns.items()}
    result["time_yr"] = result["age_Myr"] * 1e6
    fh_driver(result)
    return result


def _first_wd(data, settings, until_age_Myr):
    age = data["age_Myr"]
    if age[0] != 0 or not np.isclose(data["mass_MSun"][0], settings["initial_mass_MSun"],
                                     rtol=1e-12, atol=1e-12):
        raise ValueError("The stellar driver must start at age zero with the configured mass")
    wd = np.isin(data["stellar_type"], FH_WD_TYPES)
    first = int(np.argmax(wd)) if wd.any() else None
    if first is not None and (first == 0 or not wd[first:].all()):
        raise ValueError("Expected an uninterrupted WD phase following the progenitor")
    if until_age_Myr is None:
        if first is None or first == len(age) - 1:
            raise ValueError("Expected post-WD coverage in the full stellar track")
        end = age[first] + settings["wd_after_yr"] / 1e6
    else:
        end = until_age_Myr
    if abs(age[-1] - end) > 1e-9:
        raise ValueError(f"The stellar history ends at {age[-1]} Myr; expected {end}")
    return float(age[first]) if first is not None else None


def _checked_reference(data, reference_path):
    """Check all five stellar columns using the notebook's age and state tolerances."""
    absolute = dict(time_yr=1e-3, age_Myr=1e-9, mass_MSun=1e-12, radius_AU=1e-14,
                    stellar_type=0.)
    errors = {}
    with np.load(reference_path, allow_pickle=False) as archive:
        for key in DRIVER_FIELDS:
            actual, expected = data[key], archive[key]
            if actual.shape != expected.shape:
                raise ValueError(f"Reference {key} has shape {expected.shape}, generated {actual.shape}")
            relative = 1e-9 if key in ("radius_AU", "mass_MSun") else 0.
            matches = np.isclose(actual, expected, rtol=relative, atol=absolute[key])
            if not np.all(matches):
                i = int(np.flatnonzero(~matches)[0])
                raise ValueError(f"Stellar reference mismatch at {key}[{i}]: "
                                 f"generated {actual[i]:.16g}, expected {expected[i]:.16g}")
            errors[key] = float(np.max(np.abs(actual - expected)))
    return errors


def verify_track(path, reference_path):
    """Compare a saved driver to the machine export; return maximum absolute errors."""
    with np.load(path, allow_pickle=False) as archive:
        data = fh_driver({key: archive[key] for key in DRIVER_FIELDS})
    return _checked_reference(data, reference_path)


def _reuse(path, settings, verify_against, until_age_Myr):
    with np.load(path, allow_pickle=False) as archive:
        if "preparation_json" not in archive:
            raise ValueError(f"{path} has no preparation record. Choose a fresh output directory.")
        record = json.loads(archive["preparation_json"].item())
        if (record.get("schema") != SCHEMA or record.get("settings") != settings
                or record.get("until_age_Myr") != until_age_Myr):
            raise ValueError(f"{path} was generated with different stellar settings. Choose a fresh directory.")
        data = fh_driver({k: archive[k] for k in (*DRIVER_FIELDS, *STELLAR_STRUCTURE_FIELDS)})
        first_wd = _first_wd(data, settings, until_age_Myr)
        if first_wd != record["first_wd_age_Myr"] or data["age_Myr"][-1] != record["end_age_Myr"]:
            raise ValueError(f"Saved stellar endpoints differ from the preparation record in {path}")
        requests = archive["requested_time_Myr"]
        if (len(requests) != record["rows"] or requests.shape != data["age_Myr"].shape
                or np.any(np.diff(requests) <= 0) or not np.allclose(
                    data["age_Myr"], requests, atol=1e-9, rtol=0)):
            raise ValueError(f"Invalid saved SeBa requests in {path}")
    if verify_against:
        _checked_reference(data, verify_against)
    print(f"Reusing {len(requests):,} stellar states in {path}", flush=True)
    return path


def prepare_track(path, settings, channel="sockets", verify_against=None, *, until_age_Myr=None):
    """Generate one shared driver, or reuse a completed track with matching inputs."""
    path = Path(path).expanduser().resolve()
    settings = normalized_settings(settings)
    if until_age_Myr is not None and (not math.isfinite(until_age_Myr) or until_age_Myr <= 0):
        raise ValueError("until_age_Myr must be a finite positive age")
    reference = Path(verify_against).expanduser().resolve() if verify_against else None
    if reference is not None and (reference == path or not reference.is_file()):
        raise ValueError("The reference must be an existing driver.npz separate from the output.")
    with case_lock(path.parent / ".locks" / (path.name + ".seba.lock")):
        return _prepare_locked(path, settings, channel, reference, until_age_Myr)


def _prepare_locked(path, settings, channel, reference, until_age_Myr):
    if path.exists():
        return _reuse(path, settings, reference, until_age_Myr)
    started = perf_counter()
    with tempfile.TemporaryDirectory(prefix=".seba-", dir=path.parent) as folder:
        # AMUSE may create solver logs in the working directory.
        with chdir(folder):
            data = _generate(settings, channel, until_age_Myr)
        errors = _checked_reference(data, reference) if reference else None
        packages = {"numpy": np.__version__}
        for name in ("amuse-framework", "amuse-seba"):
            try:
                packages[name] = version(name)
            except PackageNotFoundError:
                packages[name] = "unreported"
        first_wd = _first_wd(data, settings, until_age_Myr)
        record = dict(schema=SCHEMA, settings=settings, packages=packages, channel=channel,
                      python=platform.python_version(), rows=len(data["age_Myr"]),
                      until_age_Myr=until_age_Myr,
                      first_wd_age_Myr=first_wd,
                      end_age_Myr=float(data["age_Myr"][-1]), wall_seconds=perf_counter() - started,
                      reference_max_abs_error=errors)
        pending = Path(folder) / "stellar_driver.npz"
        with pending.open("wb") as stream:
            np.savez(stream, **data, preparation_json=json.dumps(record, allow_nan=False))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(pending, path)
    print(f"Prepared {record['rows']:,} stellar states in {(perf_counter() - started) / 60:.1f} min: "
          f"{path}", flush=True)
    return path
