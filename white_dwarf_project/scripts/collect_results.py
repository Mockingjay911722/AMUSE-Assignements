"""Collect completed orbit jobs without importing AMUSE or opening large histories."""

import csv
import json
import tempfile
from itertools import groupby
from pathlib import Path


RESULT_FIELDS = (
    "suite data_level record_type case_result case a0_AU inclination_deg result outcome "
    "n1_outcome n2_outcome n2_tight_outcome incomplete_levels integrated_states "
    "saved_history_rows end_age_Myr coupling_contact_age_difference_yr "
    "solver_contact_age_difference_yr minimum_sampled_clearance_AU sampled_emax "
    "tide_variant tidal_friction tV_star_yr initial_a_out_AU reported_level final_a_in_AU "
    "final_e_in final_q_in_AU final_Q_in_AU final_a_out_AU final_e_out final_Q_out_AU "
    "initial_e_out initial_omega_out_deg stellar_friction planetary_friction tV_planet_yr "
    "distance_to_planet_band_AU"
).split()
KEY = ("suite", "case", "data_level")
# Differences are rerun minus reference. A missing measurement stays blank.
DELTA_FIELDS = {
    "end_age_Myr": ("delta_end_age_yr", 1e6),
    "minimum_sampled_clearance_AU": ("delta_minimum_sampled_clearance_AU", 1),
    "sampled_emax": ("delta_sampled_emax", 1),
    "final_a_in_AU": ("delta_final_a_in_AU", 1),
    "final_e_in": ("delta_final_e_in", 1),
    "final_a_out_AU": ("delta_final_a_out_AU", 1),
    "final_e_out": ("delta_final_e_out", 1),
    "integrated_states": ("delta_integrated_states", 1),
    "saved_history_rows": ("delta_saved_history_rows", 1),
}


def _read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path, columns, rows):
    stream = tempfile.NamedTemporaryFile(dir=path.parent, newline="", encoding="utf-8",
                                         mode="w", delete=False)
    tmp = Path(stream.name)
    try:
        with stream:
            writer = csv.DictWriter(stream, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def _result(summary):
    if summary["outcome"] == "surface_contact" and summary["within_model_domain"]:
        return "Surface contact"
    if summary["reached_requested_end"]:
        return "Reached WD endpoint" if summary.get("reached_wd_endpoint", True) else "Reached requested endpoint"
    if summary["outcome"] == "sampled_roche_crossing":
        return "Sampled Roche crossing"
    return "Model stop"


def _completed_row(folder, job, record, summary):
    meta = _read_json(folder / "job.json")
    if (meta["case"], meta["data_level"], meta["suite"]) != (
        job["case"], job["data_level"], job["suite"]
    ):
        raise ValueError(f"Job metadata and directory disagree in {folder}")
    if meta.get("pilot_age_Myr") is not None:
        raise ValueError(f"Pilot output belongs in a separate directory: {folder}")
    if (record["n_states"] != summary["n_states"] or record["outcome"] != summary["outcome"]
            or record.get("n_saved", record["n_states"]) != summary.get("n_saved", summary["n_states"])):
        raise ValueError(f"Completed run and summary disagree in {folder}")
    spec = meta["spec"]
    cfg = record["inputs"]["cfg"]
    inc = spec["inclination_deg"]
    companion = inc is not None
    suite = job["suite"]
    level = job["data_level"]
    row = dict(
        suite=suite, data_level=level, case=job["case"], record_type="Numerical summary",
        result=_result(summary), outcome=summary["outcome"], a0_AU=cfg["a_in_AU"],
        inclination_deg=inc, integrated_states=summary["n_states"],
        saved_history_rows=summary.get("n_saved", summary["n_states"]),
        end_age_Myr=summary["end_age_Myr"], sampled_emax=summary["sampled_emax"],
        minimum_sampled_clearance_AU=summary["minimum_sampled_clearance_AU"],
        tide_variant=spec.get("tide_variant", "baseline"), tidal_friction=cfg["tides"],
        stellar_friction=cfg.get("tides_star", cfg["tides"]),
        planetary_friction=cfg.get("tides_planet", cfg["tides"]),
        tV_star_yr=None if cfg.get("stellar_tide_prescription") == 1 else cfg["tV_star_yr"],
        tV_planet_yr=cfg["tV_planet_yr"], initial_a_out_AU=cfg["a_out_AU"] if companion else None,
        initial_e_out=cfg["e_out"] if companion else None,
        initial_omega_out_deg=cfg.get("omega_out_deg") if companion else None,
    )
    # Match which endpoint columns the archived main survey and follow-ups reported.
    if suite != "survey":
        row["reported_level"] = level
        if row["result"] == "Reached WD endpoint":
            a, e = summary["final_a_in_AU"], summary["final_e_in"]
            row.update(final_a_in_AU=a, final_e_in=e, final_q_in_AU=a*(1-e), final_Q_in_AU=a*(1+e))
            if companion:
                outer = summary
                if summary.get("final_a_out_AU") is None:
                    state = _read_json(folder / "checkpoint.json")["state"]
                    outer = {"final_a_out_AU": state["a_out_AU"], "final_e_out": state["e_out"]}
                a_out, e_out = outer["final_a_out_AU"], outer["final_e_out"]
                row.update(final_a_out_AU=a_out, final_e_out=e_out, final_Q_out_AU=a_out*(1+e_out))
            if suite in ("kmt-adaptive", "kmt-envelope", "kmt-eccentric") and summary["within_model_domain"]:
                row["distance_to_planet_band_AU"] = max(1.9-a*(1+e), a*(1-e)-2.3, 0.0)
    return row


def _key(row):
    return tuple(row[col] for col in KEY)


def _compare(reference, output, jobs, rows_by_key, statuses):
    with Path(reference).open(newline="", encoding="utf-8-sig") as stream:
        archived = {}
        for row in csv.DictReader(stream):
            key = _key(row)
            if key in archived:
                raise ValueError(f"Duplicate job {key} in {reference}")
            archived[key] = row
    compared = mismatches = unmatched = numerical = 0
    rows = []
    for job in jobs:
        key = _key(job)
        ref = archived.get(key)
        new = rows_by_key.get(key)
        entry = dict(zip(KEY, key))
        entry.update(status=statuses[key], reference_record_type=ref["record_type"] if ref else "",
                     reference_outcome=ref["outcome"] if ref else "", rerun_outcome=new["outcome"] if new else "")
        if ref is None:
            unmatched += 1
        elif new:
            compared += 1
            entry["outcome_matches"] = new["outcome"] == ref["outcome"]
            mismatches += int(not entry["outcome_matches"])
            if ref["record_type"] == "Numerical summary":
                numerical += 1
            for field, (col, factor) in DELTA_FIELDS.items():
                if ref.get(field, "") != "" and new.get(field) is not None:
                    entry[col] = (float(new[field]) - float(ref[field])) * factor
        rows.append(entry)
    path = output / "reference_comparison.csv"
    columns = [*KEY, "status", "reference_record_type", "reference_outcome", "rerun_outcome",
               "outcome_matches", *(col for col, factor in DELTA_FIELDS.values())]
    _write_csv(path, columns, rows)
    return dict(comparison=path, compared=compared, outcome_mismatches=mismatches,
                unmatched=unmatched, numerical_comparisons=numerical)


def collect_results(output_root, jobs, *, reference=None):
    """Write one numerical row per completed job, plus an optional independent comparison.

    Pass the configured campaign jobs to include all available levels when determining
    case agreement. Missing and interrupted jobs never become simulated measurements.
    """
    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    path = output / "planet_survival_results.csv"
    if reference is not None and Path(reference).resolve() in {
        path.resolve(), (output / "reference_comparison.csv").resolve()
    }:
        raise ValueError("Choose a separate output directory for regenerated results and the reference CSV.")
    jobs = list(jobs)
    rows_by_key, statuses = {}, {}
    for _, group in groupby(sorted(jobs, key=lambda j: (j["suite"], j["case"])),
                            key=lambda j: (j["suite"], j["case"])):
        summaries, rows, unfinished = {}, {}, []
        for job in group:
            key = _key(job)
            if key in statuses:
                raise ValueError(f"Duplicate job {key}")
            suite, case, level = key
            folder = output / "runs" / suite / case / level
            record_path = folder / "run.json"
            status = "pending"
            if record_path.is_file():
                record = _read_json(record_path)
                status = record["status"]
                if status == "completed":
                    summary = _read_json(folder / "summary.json")
                    summaries[level] = summary
                    rows[level] = _completed_row(folder, job, record, summary)
                    rows_by_key[key] = rows[level]
            elif (folder / "job.json").is_file():
                status = "initializing"
            statuses[key] = status
            if status != "completed":
                unfinished.append(f"{level}={status}")

        outcomes = {f"{level}_outcome": summary["outcome"] for level, summary in summaries.items()}
        disagree = len(set(outcomes.values())) > 1
        for level, row in rows.items():
            row.update({k: v for k, v in outcomes.items() if k in RESULT_FIELDS},
                       incomplete_levels=" ".join(unfinished),
                       case_result="Settings disagree" if disagree else row["result"])
        for first, second, field in [
            ("n1", "n2", "coupling_contact_age_difference_yr"),
            ("n2", "n2_tight", "solver_contact_age_difference_yr"),
        ]:
            if all(summaries.get(lvl, {}).get("outcome") == "surface_contact" for lvl in (first, second)):
                rows[second][field] = abs(summaries[first]["end_time_yr"] - summaries[second]["end_time_yr"])

    _write_csv(path, RESULT_FIELDS, (rows_by_key[_key(j)] for j in jobs if _key(j) in rows_by_key))
    pending = sum(s == "pending" for s in statuses.values())
    stats = dict(csv=path, completed=len(rows_by_key), pending=pending,
                 incomplete=len(jobs)-len(rows_by_key)-pending, comparison=None,
                 compared=0, outcome_mismatches=0, unmatched=0, numerical_comparisons=0)
    if reference is not None:
        stats.update(_compare(reference, output, jobs, rows_by_key, statuses))
    return stats
