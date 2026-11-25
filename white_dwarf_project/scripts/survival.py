"""Prepare, run and collect an AMUSE white-dwarf planet-survival campaign."""
import argparse
from collections import Counter
from hashlib import sha256
import json
import math
import multiprocessing
from multiprocessing.connection import wait
import os
from pathlib import Path
import signal
import sys
from time import perf_counter

from campaign import build_jobs, DEFAULT_CONFIG, load_config, select_jobs, settings_for

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT / "data" / "survival_campaign"
REFERENCE_CSV = PROJECT / "reference" / "planet_survival_results.csv"
REFERENCE_DRIVER_SHA256 = "340e5e0c992df4929dbc097d83a5a297e9ddb00899d4ff6eaeef096592f1049b"
HISTORY_ROW_BYTES = 47 * 8  # 47 float64 columns in an orbital snapshot


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("Expected a positive integer.")
    return value


def worker_count(value):
    return value if value == "auto" else positive_int(value)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("plan", "List scheduled cases and estimate storage with an optional driver."),
        ("prepare", "Generate and check the common SeBa history."),
        ("pilot", "Benchmark selected orbits on a short stellar history."),
        ("run", "Run or resume the selected case–resolution jobs."),
        ("collect", "Write the results CSV from completed production jobs."),
    ):
        cmd = sub.add_parser(name, help=help_text, description=help_text)
        cmd.add_argument("--config", type=Path, help="Campaign TOML (default: configs/survival_reference.toml).")
        cmd.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Campaign output directory.")
        if name in ("plan", "run", "pilot"):
            cmd.add_argument("--suite", nargs="+", metavar="NAME", help="Select one or more suites.")
            cmd.add_argument("--case", nargs="+", metavar="NAME", help="Select one or more exact case labels.")
            cmd.add_argument("--level", nargs="+", metavar="NAME",
                             help="Select scheduled coupling/solver levels.")
            cmd.add_argument("--history-stride", type=positive_int, metavar="N",
                             help="Save every Nth coupling state. Default comes from each selected suite.")
            cmd.add_argument("--driver", type=Path, help="Use an existing stellar driver NPZ.")
        if name in ("prepare", "pilot", "run"):
            cmd.add_argument("--channel", choices=("sockets", "mpi"), help="AMUSE channel (reference: sockets).")
        if name in ("pilot", "run"):
            cmd.add_argument("--jobs", default="auto", type=worker_count, metavar="N|auto",
                             help="Concurrent case–resolution jobs on this machine.")
            cmd.add_argument("--max-outputs", type=positive_int, default=25_000_000,
                             help="Maximum planned coupling states for one job (default: 25000000).")
        if name == "pilot":
            cmd.add_argument("--age-myr", type=float, default=1.0, help="Pilot endpoint in Myr (default: 1).")
        if name == "plan":
            cmd.add_argument("--list", action="store_true", help="Print each selected case and resolution.")
        if name == "prepare":
            cmd.add_argument("--verify-against", type=Path, metavar="DRIVER.NPZ",
                             help="Check all five stellar columns against an independent reference.")
        if name == "collect":
            cmd.add_argument("--compare", nargs="?", const=REFERENCE_CSV, type=Path, metavar="REFERENCE.CSV",
                             help="Write signed differences from the bundled or specified reference CSV.")
    return parser.parse_args()


def available_jobs():
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    quota = Path("/sys/fs/cgroup/cpu.max")
    if quota.is_file():
        cpu, period = quota.read_text().split()
        if cpu != "max":
            cpus = min(cpus, max(1, int(cpu) // int(period)))
    if os.getenv("SLURM_CPUS_PER_TASK", "").isdigit():
        cpus = min(cpus, max(1, int(os.environ["SLURM_CPUS_PER_TASK"])))
    # Limit simultaneous full drivers, target arrays and solver processes.
    memory = None
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        for line in meminfo.read_text().splitlines():
            if line.startswith("MemAvailable:"):
                memory = int(line.split()[1]) * 1024
                break
    cg = Path("/sys/fs/cgroup")
    if (cg / "memory.max").is_file() and (cg / "memory.current").is_file():
        limit = (cg / "memory.max").read_text().strip()
        if limit != "max":
            remaining = max(0, int(limit) - int((cg / "memory.current").read_text()))
            memory = remaining if memory is None else min(memory, remaining)
    memory_slots = max(1, int((memory - 1024**3) // (2 * 1024**3))) if memory is not None else 2
    return max(1, min(cpus, memory_slots))


def plan(jobs, settings, driver_path, list_jobs):
    counts_by_step = {}
    if driver_path is not None:
        from orbit_worker import load_driver
        from survey_orbits import fh_coupling_counts
        track = load_driver(driver_path)
        for step in {job["spec"]["step_yr"] for job in jobs}:
            counts, _ = fh_coupling_counts(track, step, settings.get("max_mass_fraction", 0.002),
                                           settings.get("max_radius_fraction", 0.01))
            counts_by_step[step] = int(counts.sum())
        phase_changes = int((track["stellar_type"][1:] != track["stellar_type"][:-1]).sum())

    print(f"{len(jobs)} scheduled jobs, {len({(j['suite'], j['case']) for j in jobs})} configurations")
    for suite, n in Counter(j["suite"] for j in jobs).items():
        subset = [j for j in jobs if j["suite"] == suite]
        levels = Counter(j["data_level"] for j in subset)
        print(f"  {suite:<18} {n:>3} jobs  "
              + ", ".join(f"{v} {k}" for k, v in levels.items())
              + f"  |  save every {subset[0]['history_stride']}")
    total_steps = total_rows = 0
    for j in jobs:
        extra = ""
        if counts_by_step:
            count = counts_by_step[j["spec"]["step_yr"]] * j["level"] + 1
            rows = min(count, (count-1) // j["history_stride"] + phase_changes + 2)
            total_steps += count
            total_rows += rows
            extra = f"  targets={count:,}  history<={rows * HISTORY_ROW_BYTES / 2**20:.1f} MiB"
        if list_jobs:
            print(f"{j['suite']}  {j['case']}  {j['data_level']}  stride={j['history_stride']}{extra}")
    if counts_by_step:
        print(f"All jobs reaching the supplied endpoint: {total_steps:,} coupling states; "
              f"up to {total_rows:,} history rows ({total_rows * HISTORY_ROW_BYTES / 2**30:.2f} GiB).")
    else:
        print("Add --list for case labels. Add --driver PATH for full-duration target and storage budgets.")


def generated_track(root, config, channel, *, reference=None, age=None):
    from stellar_track import prepare_track
    return prepare_track(root / "stellar_driver.npz", config["stellar"], channel=channel,
                         verify_against=reference, until_age_Myr=age)


def check_driver(path, jobs, config, age):
    """Check the input before scheduling expensive jobs, including Hurley fields."""
    import numpy as np
    from orbit_worker import load_driver
    from survey_orbits import FH_WD_TYPES
    envelope = any(j["spec"]["cfg"].get("stellar_tide_prescription") == 1 for j in jobs)
    data = load_driver(path, evolving_tide=envelope, pilot_age_Myr=age)
    if not math.isclose(data["mass_MSun"][0], config["stellar"]["initial_mass_MSun"], rel_tol=1e-12):
        raise ValueError("Initial stellar mass differs between the supplied driver and the TOML.")
    if age is None:
        first = int(np.argmax(np.isin(data["stellar_type"], FH_WD_TYPES)))
        duration = data["age_Myr"][-1] - data["age_Myr"][first]
        if not math.isclose(duration * 1e6, config["stellar"]["wd_after_yr"], abs_tol=0.001):
            raise ValueError("WD duration differs between the supplied driver and the TOML.")
    with np.load(path, allow_pickle=False) as archive:
        if "preparation_json" in archive:
            from stellar_track import normalized_settings
            prep = json.loads(archive["preparation_json"].item())
            if prep["settings"] != normalized_settings(config["stellar"]):
                raise ValueError("The supplied prepared driver's stellar settings differ from this TOML.")
            return True
    return False


def bind_campaign(root, config, driver, age, has_preparation_record):
    from orbit_worker import case_lock
    from survey_orbits import fh_json
    digest = sha256()
    with driver.open("rb") as source:
        for block in iter(lambda: source.read(2**20), b""):
            digest.update(block)
    fingerprint = digest.hexdigest()
    if fingerprint == REFERENCE_DRIVER_SHA256:
        from stellar_track import normalized_settings
        if normalized_settings(config["stellar"]) != normalized_settings({}):
            raise ValueError("This is the original machine driver. Restore the reference stellar settings "
                             "or omit --driver to evolve the new star.")
        origin = "original notebook 06 driver"
    elif has_preparation_record:
        origin = "SeBa settings verified from the preparation record"
    else:
        origin = "external driver without stellar preparation metadata"
        print("External stellar file: metallicity and sampling settings cannot be verified from its columns.")
    record = dict(schema=1, config=config, driver_sha256=fingerprint, stellar_input=origin, pilot_age_Myr=age)
    with case_lock(root / ".campaign.lock"):
        path = root / "campaign.json"
        if path.is_file():
            if json.loads(path.read_text(encoding="utf-8")) != record:
                raise ValueError(f"The campaign settings or driver changed. Choose a fresh --output directory: {root}")
        else:
            if (root / "runs").exists():
                raise ValueError(f"Runs already exist without a campaign record in {root}. Choose a fresh --output.")
            fh_json(path, record)


def execute_job(task):
    # A private POSIX group lets the parent stop this job's native AMUSE worker too.
    if os.name == "posix":
        os.setpgrp()
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    from orbit_worker import run_job
    root, driver, job, settings, max_outputs, age = task
    key = (job["suite"], job["case"], job["data_level"])
    target = root.joinpath("runs", *key)
    try:
        run_job(target, driver, job["spec"], settings, job["data_level"],
                job["history_stride"], suite=job["suite"], max_outputs=max_outputs,
                pilot_age_Myr=age)
    except Exception as exc:
        print(f"{'/'.join(key)}: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)


def stop_job(process, force=False):
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL if force else signal.SIGTERM)
            return
        except ProcessLookupError:
            pass  # The child can still be starting up, before setpgrp().
    if process.is_alive():
        if force:
            process.kill()
        else:
            process.terminate()


def schedule_jobs(tasks, workers, total):
    """Wait for child exits directly, including abrupt exits without a result."""
    context = multiprocessing.get_context("spawn")
    active = {}
    completed = failures = 0
    interrupted = 0

    def request_stop(signum, frame):
        nonlocal interrupted
        interrupted = signum

    previous = {sig: signal.signal(sig, request_stop) for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        while completed < total:
            while not interrupted and len(active) < workers:
                task = next(tasks, None)
                if task is None:
                    break
                job = task[2]
                key = (job["suite"], job["case"], job["data_level"])
                process = context.Process(target=execute_job, args=(task,))
                process.start()
                active[process.sentinel] = (process, task[0].joinpath("runs", *key), "/".join(key))
            if interrupted:
                raise SystemExit(128 + interrupted)
            # A signal delivered only to the parent also needs to wake this loop.
            for sentinel in wait(active, timeout=0.25):
                process, folder, label = active[sentinel]
                process.join()
                completed += 1
                try:
                    if process.exitcode != 0:
                        stop_job(process)  # Release any surviving native subprocess.
                        raise RuntimeError(f"worker exit code {process.exitcode}")
                    # The worker commits run.json last, after its summary and final checkpoint.
                    record = json.loads((folder / "run.json").read_text(encoding="utf-8"))
                    if record["status"] != "completed":
                        raise ValueError("worker exited before writing a completed run record")
                    print(f"[{completed}/{total}] {label}  {record['outcome']}  "
                          f"{record['n_states']:,} integrated / {record['n_saved']:,} saved", flush=True)
                except (KeyError, OSError, RuntimeError, ValueError) as exc:
                    failures += 1
                    print(f"[{completed}/{total}] FAILED {label}: {exc}", file=sys.stderr, flush=True)
                del active[sentinel]
                process.close()
    finally:
        for process, _, _ in active.values():
            stop_job(process)
        deadline = perf_counter() + 2.0
        for process, _, _ in active.values():
            process.join(timeout=max(0.0, deadline - perf_counter()))
            if process.is_alive():
                stop_job(process, force=True)
                process.join()
            process.close()
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    return failures


def simulate(args, root, config, jobs, settings):
    from importlib.util import find_spec
    if find_spec("amuse.community.secularmultiple.interface") is None:
        raise ImportError("Install the AMUSE SecularMultiple community code in this environment.")
    age = args.age_myr if args.command == "pilot" else None
    if age is not None:
        if not math.isfinite(age) or age <= 0:
            raise ValueError("--age-myr must be finite and positive.")
        root = root / "pilots" / f"{age!r}Myr"
    channel = settings["worker_options"]["channel_type"]
    driver = args.driver.expanduser().resolve() if args.driver else generated_track(root, config, channel, age=age)
    has_record = check_driver(driver, jobs, config, age)
    bind_campaign(root, config, driver, age, has_record)
    # Spawn imports the solver in each child. Numerical libraries should not
    # create a second thread pool inside each independent case.
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[name] = "1"
    workers = min(len(jobs), available_jobs() if args.jobs == "auto" else args.jobs)
    print(f"{len(jobs)} jobs, {workers} concurrent on this machine. Output: {root}", flush=True)
    started = perf_counter()
    # Finer, shorter-timestep jobs start early to reduce idle workers at the end.
    ordered = sorted(jobs, key=lambda j: j["level"]/j["spec"]["step_yr"], reverse=True)
    tasks = ((root, driver, job, settings, args.max_outputs, age) for job in ordered)
    failures = schedule_jobs(tasks, workers, len(jobs))
    print(f"Batch elapsed: {perf_counter() - started:.1f} s. {failures} failed.", flush=True)
    return 1 if failures else 0


def main():
    args = parse_args()
    root = args.output.expanduser().resolve()
    manifest_path = root / "campaign.json"
    manifest = None
    if args.command == "collect" and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("pilot_age_Myr") is not None:
            raise ValueError("Collect production results from the campaign root.")
    config = (manifest["config"] if manifest and args.config is None
              else load_config(args.config or DEFAULT_CONFIG))
    if getattr(args, "channel", None):
        config["solver"]["channel"] = args.channel
    if manifest and config != manifest["config"]:
        raise ValueError("This config differs from the one used in the output directory.")
    settings = settings_for(config)
    jobs = build_jobs(config)
    if args.command == "collect":
        from collect_results import collect_results
        stats = collect_results(root, jobs, reference=args.compare)
        print(f"{stats['completed']} completed, {stats['pending']} pending, {stats['incomplete']} incomplete. "
              f"Results: {stats['csv']}")
        if stats["comparison"]:
            print(f"{stats['compared']} outcomes compared; {stats['outcome_mismatches']} differ. "
                  f"{stats['numerical_comparisons']} numerical references. {stats['comparison']}")
        return 0
    if args.command == "prepare":
        generated_track(root, config, settings["worker_options"]["channel_type"], reference=args.verify_against)
        if args.verify_against:
            print("All five stellar columns agree with the supplied reference within the documented tolerances.")
        return 0
    if args.command == "pilot" and not any((args.suite, args.case, args.level)):
        # Keep the first run short even for a large custom config.
        jobs = [j for j in jobs if j["data_level"] == "n1"][:2] or jobs[:1]
    jobs = select_jobs(jobs, suites=args.suite, cases=args.case, levels=args.level,
                       history_stride=args.history_stride)
    if args.command == "plan":
        plan(jobs, settings, args.driver.expanduser() if args.driver else None, args.list)
        return 0
    return simulate(args, root, config, jobs, settings)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(130)
    except (ImportError, ValueError, KeyError, TypeError, RuntimeError, OSError) as exc:
        sys.exit(f"Campaign error: {exc}")
