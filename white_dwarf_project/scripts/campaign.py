"""Expand a TOML campaign into independent case–resolution jobs."""
from itertools import product
import math
from pathlib import Path
import re
import tomllib

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "survival_reference.toml"
ORBIT_G = 39.4852492465  # AU^3 / (solar mass yr^2), as in the course model


def _resolution(name):
    match = re.fullmatch(r"n([1-9][0-9]*)(_tight)?", name) if isinstance(name, str) else None
    if not match:
        raise ValueError(f"Invalid coupling level {name!r}. Use n1, n2, n2_tight, etc.")
    return int(match.group(1)), bool(match.group(2))


def load_config(path=DEFAULT_CONFIG):
    with Path(path).expanduser().open("rb") as stream:
        return tomllib.load(stream)


def _number(value, name, *, positive=False):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive.")
    return value


def _bool(value, name):
    if type(value) is not bool:
        raise ValueError(f"{name} must be true or false.")
    return value


def _stride(value):
    if type(value) is not int or value < 1:
        raise ValueError("history_stride must be a positive integer.")
    return value


def _expect_keys(data, allowed, table):
    unexpected = set(data) - set(allowed)
    if unexpected:
        raise ValueError(f"Unknown key(s) in {table}: {', '.join(sorted(unexpected))}")


def settings_for(config):
    """Return the orbital settings used by the notebook 06 integrator."""
    orbit = dict(config["orbit"])
    positive = ("m_planet_MSun", "m_bd_MSun", "a_in_AU", "a_out_AU", "r_planet_AU",
                "r_bd_AU", "kappa_star", "kappa_planet",
                "tV_star_yr", "tV_planet_yr")
    spins = ("spin_star_over_n", "spin_planet_over_n")
    angles = ("outer_inclination_deg", "mutual_inclination_deg")
    tide_flags = ("tides", "tidal_bulges", "rotation_bulges")
    _expect_keys(orbit, (*positive, *spins, *angles, *tide_flags, "k_star", "k_planet",
                         "e_in", "e_out"), "[orbit]")
    for name in positive:
        _number(orbit[name], name, positive=True)
    for name in spins:
        _number(orbit[name], name)
    for name in angles:
        if not 0 <= _number(orbit[name], name) <= 180:
            raise ValueError(f"{name} must be between 0 and 180 degrees.")
    for name in ("e_in", "e_out"):
        if not 0 < _number(orbit[name], name) < 1:
            raise ValueError(f"{name} must be between 0 and 1, excluding both endpoints.")
    for name in ("k_star", "k_planet"):
        if _number(orbit[name], name) < 0:
            raise ValueError(f"{name} must be nonnegative.")
    for name in tide_flags:
        _bool(orbit[name], name)

    physics = config["physics"]
    _expect_keys(physics, ("GR", "octupole", "wind_spin"), "[physics]")
    for name in ("GR", "octupole", "wind_spin"):
        _bool(physics[name], name)
    solver = config["solver"]
    grid_settings = {"tight_factor": 10.0, "max_mass_fraction": 0.002, "max_radius_fraction": 0.01}
    _expect_keys(solver, ("rtol", "atol", "channel", *grid_settings), "[solver]")
    numerics = {key: _number(solver.get(key, value), key, positive=True)
                for key, value in grid_settings.items()}
    if numerics["tight_factor"] <= 1:
        raise ValueError("tight_factor must exceed one to tighten solver tolerances.")
    if solver["channel"] not in ("sockets", "mpi"):
        raise ValueError("solver.channel must be sockets or mpi.")
    limits = config["limits"]
    limit_flags = ("roche_enabled", "roche_stop")
    thresholds = ("adiabatic_ratio", "spin_breakup_fraction", "double_averaging_ratio",
                  "mass_transfer_relative", "radius_transfer_relative", "radius_transfer_absolute_AU",
                  "roche_coefficient", "checkpoint_seconds", "progress_seconds")
    _expect_keys(limits, (*limit_flags, *thresholds), "[limits]")
    for name in limit_flags:
        _bool(limits[name], name)
    for name in thresholds:
        _number(limits[name], name, positive=True)
    if limits["roche_stop"] and not limits["roche_enabled"]:
        raise ValueError("roche_stop requires roche_enabled = true.")
    return dict(cfg=orbit, physics={**physics, "companion": True},
                rtol=_number(solver["rtol"], "rtol", positive=True),
                atol=_number(solver["atol"], "atol", positive=True),
                limits=dict(limits), worker_options={"channel_type": solver["channel"]}, **numerics)


def _lk_step(cfg, mass):
    """A quarter of the initial Lidov–Kozai time in years, in course code units."""
    inner = mass + cfg["m_planet_MSun"]
    total = inner + cfg["m_bd_MSun"]
    p_in = 2 * math.pi * math.sqrt(cfg["a_in_AU"]**3 / (ORBIT_G * inner))
    p_out = 2 * math.pi * math.sqrt(cfg["a_out_AU"]**3 / (ORBIT_G * total))
    lk = (8 / (15 * math.pi)) * (total / cfg["m_bd_MSun"]) * p_out**2 / p_in
    return lk * (1 - cfg["e_out"]**2)**1.5 / 4


def _tide_variants(group):
    """Return case-label suffixes and the corresponding changes to the tide model."""
    if "tV_star_yr" in group:
        if "tides" in group:
            raise ValueError("Choose either tides or tV_star_yr in each group.")
        for years in group["tV_star_yr"]:
            _number(years, "tV_star_yr", positive=True)
            yield f"star_tv{years:g}", {"tides": True, "tV_star_yr": years}
    else:
        overrides = {
            "baseline": {},
            "weak_star": {"tides": True},
            "friction_off": {"tides": False},
            "planet_tides_only": {"tides": True, "tides_star": False},
            "hurley": {"tides": True, "stellar_tide_prescription": 1, "tV_star_yr": 300.0},
        }
        for mode in group.get("tides", ["baseline"]):
            if mode not in overrides:
                raise ValueError(f"Unknown tidal mode {mode!r}. Choose from {list(overrides)}.")
            yield mode, overrides[mode]


def _group_specs(group, base, mass, base_step):
    fields = ("a_in_AU", "a_out_AU", "inclinations_deg", "e_out", "omega_out_deg",
              "tides", "tV_star_yr")
    _expect_keys(group, (*fields, "isolated"), "suite group")
    for key in fields:
        if key in group and (not isinstance(group[key], list) or not group[key]):
            raise ValueError(f"Group {key} must be a nonempty array.")
    inclinations = list(group.get("inclinations_deg", []))
    for angle in inclinations:
        if not 0 <= _number(angle, "inclinations_deg") <= 180:
            raise ValueError("Mutual inclinations must lie between 0 and 180 degrees.")
    if _bool(group.get("isolated", False), "isolated"):
        inclinations.insert(0, None)
    if not inclinations:
        raise ValueError("Set isolated = true or provide inclinations_deg.")
    variants = list(_tide_variants(group))
    for a_in, a_out, e_out, periapse, inc, (variant, override) in product(
        group.get("a_in_AU", [base["a_in_AU"]]), group.get("a_out_AU", [base["a_out_AU"]]),
        group.get("e_out", [base["e_out"]]), group.get("omega_out_deg", [0.0]),
        inclinations, variants,
    ):
        _number(a_in, "a_in_AU", positive=True)
        _number(a_out, "a_out_AU", positive=True)
        _number(periapse, "omega_out_deg")
        if not 0 < _number(e_out, "e_out") < 1:
            raise ValueError("The outer eccentricity must be between 0 and 1, excluding both endpoints.")
        cfg = {**base, "a_in_AU": a_in, "a_out_AU": a_out, "e_out": e_out, **override}
        if variant == "weak_star":
            cfg["tV_star_yr"] *= 10.0
        if "omega_out_deg" in group or "e_out" in group:
            cfg["omega_out_deg"] = periapse
        for body in ("star", "planet"):
            cfg[f"spin_{body}_over_n"] *= (a_in / base["a_in_AU"])**1.5

        # Isolated KMT controls with an explicit outer distance used the same
        # grid and initial orientation as their coplanar triples.
        step = base_step
        if inc is not None or "a_out_AU" in group:
            cfg["mutual_inclination_deg"] = 0.0 if inc is None else inc
            step = min(base_step, _lk_step(cfg, mass))

        angle = "isolated" if inc is None else f"i{inc:g}"
        label = f"a{a_in:g}_{angle}"
        if inc is not None and "a_out_AU" in group:
            label += f"_outer{a_out:g}"
        if inc is not None and ("e_out" in group or "omega_out_deg" in group):
            label += f"_e{e_out:g}_w{periapse:g}"
        if "tides" in group or "tV_star_yr" in group:
            label += f"_{variant}"
        yield dict(label=label, cfg=cfg, a0_AU=a_in, inclination_deg=inc,
                   step_yr=step, tide_variant=variant)


def build_jobs(config, *, stellar_mass=None):
    settings = settings_for(config)
    mass = config["stellar"]["initial_mass_MSun"] if stellar_mass is None else stellar_mass
    _number(mass, "initial stellar mass", positive=True)
    base = settings["cfg"]
    base_step = _lk_step(base, mass)
    jobs = []
    for suite, table in config["suites"].items():
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", suite):
            raise ValueError(f"Invalid suite name {suite!r}. Use letters, numbers, _ or -.")
        _expect_keys(table, ("levels", "refine", "history_stride", "groups"), f"[suites.{suite}]")
        stride = _stride(table["history_stride"])
        specs = [spec for group in table["groups"] for spec in _group_specs(group, base, mass, base_step)]
        labels = {spec["label"] for spec in specs}
        if len(specs) != len(labels) or not specs:
            raise ValueError(f"Suite {suite} must define a nonempty grid with unique case labels.")
        levels = table["levels"]
        refine = table.get("refine", {})
        if not isinstance(levels, list) or not levels or len(set(levels)) != len(levels):
            raise ValueError(f"{suite}.levels must be a nonempty list of unique coupling levels.")
        resolutions = {name: _resolution(name) for name in [*levels, *refine]}
        for level, names in refine.items():
            if not isinstance(names, list) or set(names) - labels:
                raise ValueError(f"Unknown case or invalid case list in {suite}.refine.{level}.")
        for spec in specs:
            for level, (subdivisions, tight) in resolutions.items():
                if level in levels or spec["label"] in refine.get(level, []):
                    jobs.append(dict(suite=suite, case=spec["label"], data_level=level,
                                     level=subdivisions, tight=tight, history_stride=stride, spec=spec))
    if not jobs:
        raise ValueError("The campaign must include at least one case.")
    return jobs


def select_jobs(jobs, *, suites=None, cases=None, levels=None, history_stride=None):
    """Filter the scheduled jobs. A stride override only changes the saved histories."""
    selected = list(jobs)
    for field, values in (("suite", suites), ("case", cases), ("data_level", levels)):
        if values:
            wanted = {values} if isinstance(values, str) else set(values)
            missing = wanted - {job[field] for job in selected}
            if missing:
                raise ValueError(f"Unknown or unscheduled {field} in selection: {', '.join(sorted(missing))}")
            selected = [job for job in selected if job[field] in wanted]
    if history_stride is not None:
        stride = _stride(history_stride)
        selected = [{**job, "history_stride": stride} for job in selected]
    return selected
