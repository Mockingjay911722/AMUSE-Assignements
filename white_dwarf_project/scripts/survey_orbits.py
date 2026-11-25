"""Orbital integration and checked restarts from notebook 06."""
import hashlib
import json
import os
import warnings
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

import numpy as np

from campaign import ORBIT_G as G

STELLAR_STRUCTURE_FIELDS = (
    "luminosity_LSun", "convective_envelope_mass_MSun", "convective_envelope_radius_RSun",
)


@contextmanager
def coupled_working_directory(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


class TidalAmuseBackend:
    """Keep AMUSE calls separate from the interval and comparison logic."""
    def __init__(self, cfg, mass, radius, quadrupole, rtol, atol, events, worker_options):
        from amuse.units import units
        from amuse.datamodel import Particles
        from amuse.community.secularmultiple.interface import SecularMultiple
        self.units = units
        self.cfg = cfg
        self.solver = SecularMultiple(**worker_options)
        self.H = 0.0
        try:
            p = self.solver.parameters
            p.relative_tolerance = rtol
            p.absolute_tolerance_eccentricity_vectors = atol
            p.include_quadrupole_order_terms = quadrupole
            for name in ("include_octupole_order_binary_pair_terms", "include_octupole_order_binary_triplet_terms",
                         "include_hexadecupole_order_binary_pair_terms", "include_dotriacontupole_order_binary_pair_terms"):
                setattr(p, name, False)
            tree = Particles(5)
            tree.is_binary = [False, False, False, True, True]
            tree.child1 = None
            tree.child2 = None
            tree.mass = [mass, cfg["m_planet_MSun"], cfg["m_bd_MSun"], 0, 0] | units.MSun
            tree.radius = [radius, cfg["r_planet_AU"], cfg["r_bd_AU"], 1e-15, 1e-15] | units.au
            for axis in "xyz":
                setattr(tree, "spin_vec_"+axis, 0.0 | (1/units.yr))
            tree.semimajor_axis = 1.0 | units.au
            tree.eccentricity = 1e-3
            tree.inclination = 0.1
            tree.argument_of_pericenter = 0.0
            tree.longitude_of_ascending_node = 0.0
            tree[3].child1, tree[3].child2 = tree[0], tree[1]
            tree[4].child1, tree[4].child2 = tree[3], tree[2]
            for suffix, i, inc in (("in", 3, cfg["outer_inclination_deg"]+cfg["mutual_inclination_deg"]),
                                   ("out", 4, cfg["outer_inclination_deg"])):
                tree[i].semimajor_axis = cfg["a_"+suffix+"_AU"] | units.au
                tree[i].eccentricity = cfg["e_"+suffix]
                tree[i].inclination = float(np.deg2rad(inc))
            tree[4].argument_of_pericenter = float(np.deg2rad(cfg.get("omega_out_deg", 0.0)))
            for name in ("include_pairwise_1PN_terms", "include_pairwise_25PN_terms", "include_tidal_friction_terms",
                         "include_tidal_bulges_precession_terms", "include_rotation_precession_terms",
                         "check_for_secular_breakdown", "check_for_dynamical_instability",
                         "check_for_physical_collision_or_orbit_crossing", "check_for_minimum_periapse_distance",
                         "check_for_RLOF_at_pericentre", "sample_orbital_phases_randomly"):
                setattr(tree, name, False)
            tree.mass_dot_external = 0.0 | (units.MSun/units.yr)
            tree.radius_dot_external = 0.0 | (units.au/units.yr)
            tree.radius_ddot_external = 0.0 | (units.au/units.yr**2)
            if events:
                tree[3].check_for_physical_collision_or_orbit_crossing = True
                tree[4].check_for_physical_collision_or_orbit_crossing = True
                tree[3].check_for_secular_breakdown = quadrupole
                tree[3].check_for_dynamical_instability = True
                tree[3].dynamical_instability_criterion = 0
            normal = np.array([0., -np.sin(tree[3].inclination), np.cos(tree[3].inclination)])
            frequency = np.sqrt(G*(mass+cfg["m_planet_MSun"])/cfg["a_in_AU"]**3)
            for i, name in ((0, "star"), (1, "planet")):
                tree[i].tides_method = 2
                tree[i].tides_viscous_time_scale_prescription = 0
                tree[i].tides_viscous_time_scale = cfg[f"tV_{name}_yr"] | units.yr
                tree[i].tides_apsidal_motion_constant = cfg[f"k_{name}"]
                tree[i].tides_gyration_radius = cfg[f"kappa_{name}"]
                tree[i].include_tidal_friction_terms = cfg.get(f"tides_{name}", cfg["tides"])
                tree[i].include_tidal_bulges_precession_terms = cfg["tidal_bulges"]
                tree[i].include_rotation_precession_terms = cfg["rotation_bulges"]
                for axis, value in zip("xyz", normal*frequency*cfg[f"spin_{name}_over_n"]):
                    setattr(tree[i], "spin_vec_"+axis, float(value) | (1/units.yr))
            self.remote = self.solver.particles.add_particles(tree)
            tree.new_channel_to(self.solver.particles).copy()
            self.solver.commit_particles()
            self.H = float(self.solver.initial_hamiltonian.value_in(units.MSun*units.au**2/units.yr**2))
        except BaseException:
            self.stop()
            raise

    def set_rates(self, mdot, rdot):
        self.remote[0].mass_dot_external = float(mdot) | (self.units.MSun/self.units.yr)
        self.remote[0].radius_dot_external = float(rdot) | (self.units.au/self.units.yr)

    def advance(self, dt):
        u = self.units
        elapsed, H, flag, error = self.solver.evolve_interface(0.0 | u.yr, float(dt) | u.yr)
        elapsed = float(elapsed.value_in(u.yr))
        self.H = float(H.value_in(u.MSun*u.au**2/u.yr**2))
        flag, error = int(flag), int(error)
        if error != 0 or flag not in (0, 2):
            raise RuntimeError(f"SecularMultiple flag={flag}, error={error}")
        events = []
        if flag == 2:
            for i, name in ((3, "surface_contact"), (4, "orbit_crossing")):
                if self.remote[i].physical_collision_or_orbit_crossing_has_occurred:
                    events.append(name)
            for attribute, name in (("secular_breakdown_has_occurred", "secular_breakdown"),
                                    ("dynamical_instability_has_occurred", "instability_criterion")):
                if getattr(self.remote[3], attribute):
                    events.append(name)
            if not events:
                raise RuntimeError("Root return without an identified event.")
        return elapsed, "+".join(events)

    def snapshot(self):
        u, p = self.units, self.remote
        row = {"mass_MSun": float(p[0].mass.value_in(u.MSun)),
               "radius_AU": float(p[0].radius.value_in(u.au)), "hamiltonian": self.H}
        normals = []
        for suffix, i in (("in", 3), ("out", 4)):
            row["a_"+suffix+"_AU"] = float(p[i].semimajor_axis.value_in(u.au))
            row["e_"+suffix] = float(p[i].eccentricity)
            h = np.asarray(p[i].h_vec.value_in(u.MSun*u.au**2/u.yr), dtype=float)
            for axis, value in zip("xyz", h):
                row[f"h_{suffix}_{axis}"] = float(value)
            normals.append(h/np.linalg.norm(h))
        row["mutual_inclination_deg"] = float(np.rad2deg(np.arccos(np.clip(np.dot(*normals), -1, 1))))
        for i, name in ((0, "star"), (1, "planet")):
            for axis in "xyz":
                row[f"omega_{name}_{axis}"] = float(
                    getattr(p[i], "spin_vec_"+axis).value_in(1/u.yr)
                )
        return row

    def stop(self):
        try:
            self.solver.stop()
        except Exception as exc:
            warnings.warn(f"Worker cleanup: {exc!r}")


class WDCheckpointBackend(TidalAmuseBackend):
    def checkpoint(self):
        state = self.snapshot()
        for suffix, i in (('in', 3), ('out', 4)):
            for axis in 'xyz':
                state[f'ev_{suffix}_{axis}'] = float(getattr(self.remote[i], 'e_vec_'+axis))
        return state

    def restore(self, state):
        u = self.units
        particles = self.solver.particles
        self.remote[0].mass = state['mass_MSun'] | u.MSun
        self.remote[0].radius = state['radius_AU'] | u.au

        # AMUSE requires all six orbital components in one setter call.
        names = [f'{prefix}_vec_{axis}' for prefix in ('e', 'h') for axis in 'xyz']
        values = [np.array([state[f'ev_{suffix}_{axis}'] for suffix in ('in', 'out')])
                  for axis in 'xyz']
        values += [np.array([state[f'h_{suffix}_{axis}'] for suffix in ('in', 'out')])
                   | (u.MSun*u.au**2/u.yr) for axis in 'xyz']
        indices = particles.get_indices_of_keys(self.remote[[3, 4]].key)
        particles.set_values_in_store(indices, names, values)

        # The spin setter likewise expects all three components together.
        names = [f'spin_vec_{axis}' for axis in 'xyz']
        values = [np.array([state[f'omega_{name}_{axis}'] for name in ('star', 'planet')])
                  | (1/u.yr) for axis in 'xyz']
        indices = particles.get_indices_of_keys(self.remote[[0, 1]].key)
        particles.set_values_in_store(indices, names, values)

        self.set_rates(0., 0.)
        elapsed, event = self.advance(0.)
        if event or elapsed != 0:
            raise RuntimeError('Unexpected event while restoring the checkpoint.')


def wd_state_errors(left, right, n0, full=False):
    errors = {}
    for k in ('mass_MSun', 'radius_AU', 'a_in_AU', 'a_out_AU'):
        errors[k] = float(abs(left[k]-right[k])/max(abs(right[k]), 1e-30))
    for suffix in ('in', 'out'):
        errors['e_'+suffix] = float(abs(left['e_'+suffix]-right['e_'+suffix]))
        a = np.array([left[f'h_{suffix}_{c}'] for c in 'xyz'])
        b = np.array([right[f'h_{suffix}_{c}'] for c in 'xyz'])
        errors['h_'+suffix] = float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30))
        if full:
            errors['ev_'+suffix] = float(np.linalg.norm(
                [left[f'ev_{suffix}_{c}']-right[f'ev_{suffix}_{c}'] for c in 'xyz']))
    for name in ('star', 'planet'):
        a = np.array([left[f'omega_{name}_{c}'] for c in 'xyz'])
        b = np.array([right[f'omega_{name}_{c}'] for c in 'xyz'])
        errors['spin_'+name] = float(np.linalg.norm(a-b)/max(n0, np.linalg.norm(b)))
    if not np.isfinite(list(errors.values())).all():
        raise RuntimeError('Non-finite checkpoint comparison.')
    return errors


def wd_require_state(left, right, n0, full=False):
    errors = wd_state_errors(left, right, n0, full)
    # Tighter than our scientific comparison limits; these check reproduction.
    limits = {k: (1e-8 if k.startswith(('e_', 'ev_')) else 1e-6) for k in errors}
    if any(errors[k] > limits[k] for k in errors):
        raise RuntimeError(f'Replay/restart did not reproduce the reference: {errors}')
    return errors


class TidalCheckBackend(WDCheckpointBackend):
    def __init__(self, *args, extra_physics=None, **kwargs):
        self.extra_physics = dict(extra_physics or {})
        super().__init__(*args, **kwargs)
        try:
            self.solver.parameters.include_octupole_order_binary_pair_terms = bool(
                self.extra_physics.get('octupole', False))
            self.remote[3].include_pairwise_1PN_terms = bool(
                self.extra_physics.get('GR', False))
            self.set_rates(0., 0.)
            self.advance(0.)
        except BaseException:
            self.stop()
            raise

    def set_rates(self, mdot, rdot):
        super().set_rates(mdot, rdot)
        if not self.extra_physics.get('wind_spin', False) and getattr(self,'_wind_zeroed',False):
            return
        rate = np.zeros(3)
        if self.extra_physics.get('wind_spin', False):
            if mdot > 0:
                raise ValueError('This wind prescription only covers mass loss.')
            mass = float(self.remote[0].mass.value_in(self.units.MSun))
            omega = np.array([float(getattr(self.remote[0], 'spin_vec_'+a).value_in(
                1/self.units.yr)) for a in 'xyz'])
            # Only add wind torque / I; the worker already changes I with M and R.
            rate = (2/(3*self.cfg['kappa_star']))*(mdot/mass)*omega
        particles = self.solver.particles
        indices = particles.get_indices_of_keys(self.remote[[0]].key)
        names = [f'spin_vec_{a}_dot_external' for a in 'xyz']
        values = [np.array([v]) | (1/self.units.yr**2) for v in rate]
        particles.set_values_in_store(indices, names, values)
        self._wind_zeroed = not self.extra_physics.get('wind_spin',False)


def enrich_state(row, cfg, initial_mass):
    mp = cfg['m_planet_MSun']
    spin_J = np.zeros(3)
    spin_energy = 0.0
    for name, mass, radius in (
        ('star', row['mass_MSun'], row['radius_AU']),
        ('planet', mp, cfg['r_planet_AU']),
    ):
        omega = np.array([row[f'omega_{name}_{axis}'] for axis in 'xyz'])
        inertia = cfg[f'kappa_{name}']*mass*radius**2
        spin_J += inertia*omega
        spin_energy += 0.5*inertia*np.dot(omega, omega)
        row[f'spin_{name}_over_breakup'] = float(np.linalg.norm(omega)/np.sqrt(G*mass/radius**3))
    hin = np.array([row[f'h_in_{axis}'] for axis in 'xyz'])
    hout = np.array([row[f'h_out_{axis}'] for axis in 'xyz'])
    for axis, inner, total in zip('xyz', hin+spin_J, hin+hout+spin_J):
        row[f'J_inner_spin_{axis}'] = float(inner)
        row[f'J_total_{axis}'] = float(total)
    row['E_mechanical'] = -G*row['mass_MSun']*mp/(2*row['a_in_AU'])+spin_energy
    row['q_in_AU'] = row['a_in_AU']*(1-row['e_in'])
    row['clearance_AU'] = row['q_in_AU']-row['radius_AU']-cfg['r_planet_AU']
    for suffix in ('in', 'out'):
        extra = mp+(cfg['m_bd_MSun'] if suffix == 'out' else 0)
        row[f'aM_{suffix}_ratio'] = (
            row[f'a_{suffix}_AU']*(row['mass_MSun']+extra)
            /(cfg[f'a_{suffix}_AU']*(initial_mass+extra))
        )
    return row


class FullHistoryBackend(TidalCheckBackend):
    def __init__(self, cfg, mass, radius, quadrupole, rtol, atol, events,
                 worker_options, extra_physics=None):
        physics = dict(extra_physics or {})
        companion = bool(physics.get("companion", True))
        physics["octupole"] = companion and bool(physics.get("octupole", False))
        super().__init__(cfg, mass, radius, companion and quadrupole, rtol, atol,
                         False, worker_options, extra_physics=physics)
        try:
            # A dummy outer orbit is bookkeeping only in the isolated control.
            self.solver.parameters.include_octupole_order_binary_triplet_terms = False
            self.remote[3].check_for_physical_collision_or_orbit_crossing = bool(events)
            self.remote[3].dynamical_instability_criterion = 0
            self.remote[3].check_for_dynamical_instability = events and companion
            self.remote[3].check_for_secular_breakdown = events and companion
            self.remote[4].check_for_physical_collision_or_orbit_crossing = events and companion
            elapsed, event = self.advance(0.)
            if elapsed != 0 or event:
                raise RuntimeError("Unexpected event while applying the physics switches.")
        except BaseException:
            self.stop()
            raise

    def set_stellar_structure(self, stellar_type, luminosity, envelope_mass, envelope_radius):
        u = self.units
        star = self.remote[0]
        star.stellar_type = int(stellar_type) | u.stellar_type
        star.luminosity = float(luminosity) | u.LSun
        star.convective_envelope_mass = float(envelope_mass) | u.MSun
        star.convective_envelope_radius = (float(envelope_radius)*FH_RSUN_AU) | u.au


# Stellar driver and event diagnostics

def fh_json(path, value):
    """Replace small records only after the complete new record has been written."""
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def fh_driver(driver):
    keys = ("time_yr", "mass_MSun", "radius_AU", "stellar_type", "age_Myr")
    data = {k: np.asarray(driver[k], dtype=float) for k in keys}
    t = data["time_yr"]
    if (t.ndim != 1 or len(t) < 2 or t[0] != 0
            or any(v.shape != t.shape or not np.isfinite(v).all() for v in data.values())
            or np.any(np.diff(t) <= 0) or np.any(np.diff(data["age_Myr"]) <= 0)
            or np.any(data["mass_MSun"] <= 0) or np.any(data["radius_AU"] <= 0)
            or np.any(np.diff(data["mass_MSun"]) > 1e-12)
            or not np.equal(data["stellar_type"], np.floor(data["stellar_type"])).all()):
        raise ValueError("Require finite, increasing ages, positive mass/radius, and a mass-loss track.")
    if not np.allclose(t / 1e6 + data["age_Myr"][0], data["age_Myr"], rtol=0, atol=2e-9):
        raise ValueError("Local time and stellar age disagree.")
    if any(k in driver for k in STELLAR_STRUCTURE_FIELDS):
        if not all(k in driver for k in STELLAR_STRUCTURE_FIELDS):
            raise ValueError("The stellar structure columns must be supplied together.")
        for key in STELLAR_STRUCTURE_FIELDS:
            values = np.asarray(driver[key], dtype=float)
            if values.shape != t.shape or not np.isfinite(values).all() or np.any(values < 0):
                raise ValueError(f"Invalid stellar structure column: {key}.")
            data[key] = values
        if np.any(data["luminosity_LSun"] <= 0):
            raise ValueError("Stellar luminosity must be positive for the Hurley tide.")
    return data


def fh_limits(overrides=None):
    limits = {"adiabatic_ratio": 0.1, "spin_breakup_fraction": 0.1,
              "double_averaging_ratio": 0.01,
              "mass_transfer_relative": 1e-9, "radius_transfer_relative": 1e-8,
              "radius_transfer_absolute_AU": 1e-11,
              "roche_enabled": False, "roche_stop": False, "roche_coefficient": 2.44,
              "checkpoint_seconds": 60., "progress_seconds": 30.}
    unknown = set(overrides or {}) - set(limits)
    if unknown:
        raise ValueError(f"Unknown full-history limits: {sorted(unknown)}")
    limits.update(overrides or {})
    if any(not np.isfinite(v) or v <= 0 for k, v in limits.items()
           if k not in ("roche_enabled", "roche_stop")):
        raise ValueError("All numerical limits must be positive and finite.")
    if limits["roche_stop"] and not limits["roche_enabled"]:
        raise ValueError("Enable the sampled Roche diagnostic before requesting its stop.")
    return limits


def fh_row(state, actual, driver, cfg, physics, mdot, limits):
    row = dict(state)
    t = driver["time_yr"]
    k = min(int(np.searchsorted(t, actual, side="right")-1), len(t)-1)
    row.update(time_yr=float(actual), age_Myr=float(driver["age_Myr"][0] + actual/1e6),
               driver_mass_MSun=float(np.interp(actual, t, driver["mass_MSun"])),
               driver_radius_AU=float(np.interp(actual, t, driver["radius_AU"])),
               stellar_type=int(driver["stellar_type"][k]))
    if not np.isfinite(list(row.values())).all():
        raise RuntimeError("The worker returned a non-finite state.")
    for suffix in ("in", "out"):
        if row[f"a_{suffix}_AU"] <= 0 or not 0 <= row[f"e_{suffix}"] < 1:
            raise RuntimeError(f"The worker returned an unbound or invalid {suffix} orbit.")
    row = enrich_state(row, cfg, driver["mass_MSun"][0])
    for suffix in ("in", "out"):
        mass = row["mass_MSun"] + cfg["m_planet_MSun"]
        if suffix == "out":
            mass += cfg["m_bd_MSun"]
        period = 2*np.pi*np.sqrt(row[f"a_{suffix}_AU"]**3/(G*mass))
        row[f"period_{suffix}_yr"] = float(period)
        row[f"adiabatic_ratio_{suffix}"] = float(period*abs(mdot)/mass)
    inner_mass = row["mass_MSun"] + cfg["m_planet_MSun"]
    total_mass = inner_mass + cfg["m_bd_MSun"]
    lk = ((8/(15*np.pi))*total_mass/cfg["m_bd_MSun"]
          * row["period_out_yr"]**2/row["period_in_yr"]
          * (1-row["e_out"]**2)**1.5)
    row["LK_time_yr"] = float(lk)
    row["double_averaging_ratio"] = float(row["period_out_yr"]/lk)
    row["hierarchy_ratio"] = float(row["a_out_AU"]*(1-row["e_out"])
                                   /(row["a_in_AU"]*(1+row["e_in"])))
    # An illustrative fluid Roche scale, not a disruption calculation.
    row["roche_distance_AU"] = float(limits["roche_coefficient"] * cfg["r_planet_AU"]
                                    * (row["mass_MSun"]/cfg["m_planet_MSun"])**(1/3))
    row["roche_clearance_AU"] = row["q_in_AU"] - row["roche_distance_AU"]
    return row


def fh_sampled_stop(row, cfg, physics, limits):
    if row["clearance_AU"] <= 0:
        return "surface_contact"
    if physics.get("companion", True):
        qout = row["a_out_AU"] * (1-row["e_out"])
        if qout <= row["a_in_AU"]*(1+row["e_in"]) + cfg["r_bd_AU"]:
            return "orbit_crossing"
        # Same criterion-0 expression used by the worker, also checked at t=0.
        ratio = cfg["m_bd_MSun"] / (row["mass_MSun"]+cfg["m_planet_MSun"])
        critical = (row["a_in_AU"]*2.8*((1+ratio)*(1+row["e_out"])
                    / np.sqrt(1-row["e_out"]))**0.4
                    * (1-0.3*row["mutual_inclination_deg"]/180))
        if qout <= critical:
            return "dynamical_instability"
        if row["double_averaging_ratio"] > limits["double_averaging_ratio"]:
            return "sampled_double_averaging_limit"
    active = ("in", "out") if physics.get("companion", True) else ("in",)
    if any(row[f"adiabatic_ratio_{s}"] > limits["adiabatic_ratio"] for s in active):
        return "sampled_adiabatic_model_limit"
    if any(row[f"spin_{s}_over_breakup"] > limits["spin_breakup_fraction"]
           for s in ("star", "planet")):
        return "sampled_spin_model_limit"
    if limits["roche_enabled"] and limits["roche_stop"] and row["roche_clearance_AU"] <= 0:
        return "sampled_roche_crossing"
    return ""


def fh_read_case(folder):
    folder = Path(folder)
    record = json.loads((folder/"run.json").read_text())
    if record.get("status") != "completed":
        raise RuntimeError(f"The case is incomplete: {folder}")
    table = np.load(folder/"history.npy", mmap_mode="r", allow_pickle=False)
    size = record.get("n_saved", record["n_states"])
    if table.dtype.names is None or not 1 <= size <= len(table):
        raise ValueError("Invalid saved history layout.")
    history = {k: table[k][:size] for k in table.dtype.names}
    times = history["time_yr"]
    if times[0] != 0 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Saved history must start at zero with finite, strictly increasing times.")
    summary = json.loads((folder/"summary.json").read_text())
    if (summary.get("n_saved", summary["n_states"]) != size
            or summary["n_states"] != record["n_states"]
            or summary["outcome"] != record["outcome"] or summary["end_time_yr"] != times[-1]):
        raise ValueError("Saved summary and history disagree.")
    return history, summary


# Coupling loop, output histories, and checked restarts

def fh_history_capacity(driver, target_count, stride):
    if stride == 1:
        return target_count
    changes = int(np.count_nonzero(np.diff(driver["stellar_type"])))
    # One row per stride, plus any off-stride phase changes and the endpoint.
    return min(target_count, (target_count - 1) // stride + changes + 2)


def fh_array_digest(arrays):
    digest = hashlib.sha256()
    for array in arrays:
        digest.update(str((array.shape, array.dtype.str)).encode("ascii"))
        digest.update(memoryview(np.ascontiguousarray(array)).cast("B"))
    return digest.hexdigest()


def fh_update_extrema(extrema, row, companion):
    mass = row["driver_mass_MSun"]
    radius = row["driver_radius_AU"]
    peaks = {
        "sampled_emax": row["e_in"],
        "max_rel_mass_driver_error": abs(row["mass_MSun"] / mass - 1),
        "max_rel_radius_driver_error": abs(row["radius_AU"] / radius - 1),
        "max_abs_radius_driver_error_AU": abs(row["radius_AU"] - radius),
        "max_aM_in_error": abs(row["aM_in_ratio"] - 1),
        "max_aM_out_error": abs(row["aM_out_ratio"] - 1),
        "max_spin_star_over_breakup": row["spin_star_over_breakup"],
        "max_spin_planet_over_breakup": row["spin_planet_over_breakup"],
        "max_adiabatic_ratio_in": row["adiabatic_ratio_in"],
    }
    if companion:
        peaks["max_adiabatic_ratio_out"] = row["adiabatic_ratio_out"]
        peaks["max_sampled_double_averaging_ratio"] = row["double_averaging_ratio"]
        key = "minimum_sampled_hierarchy_ratio"
        extrema[key] = min(extrema.get(key, float("inf")), row["hierarchy_ratio"])
    for key, value in peaks.items():
        extrema[key] = max(extrema.get(key, 0.0), value)
    if row["clearance_AU"] < extrema.get("minimum_sampled_clearance_AU", float("inf")):
        extrema["minimum_sampled_clearance_AU"] = row["clearance_AU"]
        extrema["age_at_minimum_clearance_Myr"] = row["age_Myr"]


def fh_run_case(folder, driver, cfg, targets, physics, rtol=1e-14, atol=1e-16,
                worker_options=None, limits=None, resume=False, backend_factory=None,
                history_stride=1, pilot=False):
    folder = Path(folder).resolve()
    driver = fh_driver(driver)
    evolving_tide = cfg.get("stellar_tide_prescription", 0) == 1
    if evolving_tide and not all(key in driver for key in STELLAR_STRUCTURE_FIELDS):
        raise ValueError("The envelope-dependent stellar tide needs the checked SeBa structure.")
    targets = np.asarray(targets, dtype=float)
    physics = {"companion": True, "GR": True, "octupole": True, "wind_spin": True, **physics}
    if not physics["companion"]:
        physics["octupole"] = False
    worker_options = dict(worker_options or {})
    limits = fh_limits(limits)
    if not isinstance(history_stride, int) or history_stride < 1:
        raise ValueError("Use an integer history stride of at least one.")
    t = driver["time_yr"]
    if (targets.ndim != 1 or len(targets) < 2 or not np.isfinite(targets).all()
            or targets[0] != 0 or targets[-1] != t[-1] or np.any(np.diff(targets) <= 0)
            or not np.array_equal(targets[np.searchsorted(targets, t)], t)):
        raise ValueError("The grid must contain every driver knot and exactly cover the driver.")
    capacity = fh_history_capacity(driver, len(targets), history_stride)
    required_positive = ("m_planet_MSun", "m_bd_MSun", "r_planet_AU", "r_bd_AU",
                         "a_in_AU", "a_out_AU", "kappa_star", "kappa_planet",
                         "tV_star_yr", "tV_planet_yr")
    if (not all(np.isfinite(v) for v in cfg.values())
            or any(cfg[k] <= 0 for k in required_positive)
            or any(not 0 < cfg[k] < 1 for k in ("e_in", "e_out"))
            or min(cfg["k_star"], cfg["k_planet"]) < 0
            or not np.isfinite([rtol, atol]).all() or min(rtol, atol) <= 0):
        raise ValueError("Invalid orbital, tidal or solver inputs.")
    environment = {}
    for package in ("amuse-framework", "amuse-secularmultiple", "numpy"):
        try:
            environment[package] = version(package)
        except Exception:
            environment[package] = "unavailable"
    inputs = {"schema": "planet_survival_orbit_v1", "cfg": cfg, "physics": physics,
              "rtol": rtol, "atol": atol, "worker_options": worker_options,
              "limits": limits, "environment": environment, "planned_outputs": len(targets)}
    inputs.update(pilot=pilot, history_stride=history_stride, driver_digest=fh_array_digest(driver.values()),
                  targets_digest=fh_array_digest((targets,)))
    old = None
    if (folder/"run.json").is_file():
        old = json.loads((folder/"run.json").read_text())
        if old.get("inputs") != inputs:
            raise RuntimeError("The existing case has different inputs, history settings or package versions.")
        if old.get("status") == "completed":
            return fh_read_case(folder)
        if not resume:
            raise RuntimeError("The case is incomplete. Keep it, then set resume=True for a checked restart.")
        if not (folder/"checkpoint.json").is_file():
            raise RuntimeError("No complete checkpoint is available; use a new case folder.")
    else:
        if (folder/"history.npy").exists() or (folder/"checkpoint.json").exists():
            raise RuntimeError(f"Found an existing trajectory without a run record in {folder}.")
        folder.mkdir(parents=True, exist_ok=True)
        fh_json(folder/"run.json", {"status": "started", "inputs": inputs, "n_states": 0})
    factory = FullHistoryBackend if backend_factory is None else backend_factory
    rates = {k: np.diff(driver[k])/np.diff(t) for k in ("mass_MSun", "radius_AU")}
    n0 = np.sqrt(G*(driver["mass_MSun"][0]+cfg["m_planet_MSun"])/cfg["a_in_AU"]**3)
    table, backend = None, None
    size, saved_size, index, actual = 0, 0, 0, 0.
    extrema = {}
    outcome, event_origin = "completed_window", "end_of_driver"
    wall_before = float((old or {}).get("wall_seconds", 0.))
    started = last_print = last_checkpoint = perf_counter()
    phase_records, roche_record, anchor = [], None, None

    def interval_rates(k):
        return float(rates["mass_MSun"][k]), float(rates["radius_AU"][k])

    def interval_index(start):
        return min(int(np.searchsorted(t, start, side="right")-1), len(t)-2)

    tide_knot = None

    def apply_stellar_tide(k):
        nonlocal tide_knot
        if k == tide_knot:
            return
        backend.set_stellar_structure(driver["stellar_type"][k], *(
            driver[key][k] for key in STELLAR_STRUCTURE_FIELDS
        ))
        tide_knot = k

    def advance_one(start, end):
        k = interval_index(start)
        if evolving_tide:
            apply_stellar_tide(k)
        mdot, rdot = interval_rates(k)
        backend.set_rates(mdot, rdot)
        elapsed, event = backend.advance(float(end-start))
        if (not np.isfinite(elapsed) or elapsed < 0 or elapsed > end-start+1e-7
                or (not event and not np.isclose(elapsed, end-start, rtol=1e-12, atol=1e-7))):
            raise RuntimeError("The worker did not reach the requested time correctly.")
        return (float(start+elapsed) if event else float(end)), event

    def checkpoint():
        nonlocal last_checkpoint
        table.flush()
        record = {"index": index, "n_states": size, "time_yr": actual,
                  "state": backend.checkpoint(), "anchor": anchor,
                  "phase_records": phase_records, "roche_record": roche_record,
                  "outcome": outcome, "event_origin": event_origin,
                  "terminal": outcome != "completed_window" or index == len(targets)-1,
                  "wall_seconds": wall_before+perf_counter()-started}
        record.update(n_saved=saved_size, row={k: float(v) for k, v in row.items()},
                      extrema=extrema, last_saved_time_yr=float(table["time_yr"][saved_size - 1]))
        fh_json(folder/"checkpoint.json", record)
        last_checkpoint = perf_counter()
        return record

    try:
        with coupled_working_directory(folder):
            if old is not None:
                cp = json.loads((folder/"checkpoint.json").read_text())
                table = np.lib.format.open_memmap(folder/"history.npy", mode="r+")
                index, size, actual = int(cp["index"]), int(cp["n_states"]), float(cp["time_yr"])
                saved_size = int(cp.get("n_saved", size))
                if (size != index+1 or not 1 <= saved_size <= size <= len(targets)
                        or len(table) != capacity or saved_size > capacity):
                    raise ValueError("Checkpoint and flushed history disagree.")
                row, extrema = cp["row"], cp["extrema"]
                times = table["time_yr"][:saved_size]
                if (row["time_yr"] != actual or times[-1] != cp["last_saved_time_yr"]
                        or times[-1] > actual or np.any(np.diff(times) <= 0)
                        or (cp["terminal"] and times[-1] != actual)
                        or (not cp["terminal"] and targets[index] != actual)):
                    raise ValueError("Checkpoint and compact history disagree.")
                scheduled_times = times[:-1] if cp["terminal"] else times
                locations = np.searchsorted(targets, scheduled_times)
                if np.any(locations >= len(targets)):
                    raise ValueError("Compact history contains out-of-range times.")
                np.testing.assert_array_equal(targets[locations], scheduled_times)
                if times[-1] == actual:
                    stored = {k: float(table[k][saved_size-1]) for k in table.dtype.names}
                    if stored != row:
                        raise ValueError("Checkpoint row and flushed history disagree.")
                wd_require_state(row, cp["state"], n0, full=False)
                phase_records, roche_record = cp["phase_records"], cp["roche_record"]
                outcome, event_origin, anchor = cp["outcome"], cp["event_origin"], cp["anchor"]
                wall_before = float(cp["wall_seconds"])
                begin = len(targets) if cp["terminal"] else index+1
            else:
                begin = 0

            if begin < len(targets):
                backend = factory(cfg, driver["mass_MSun"][0], driver["radius_AU"][0],
                                  bool(physics["companion"]), rtol, atol, True, worker_options,
                                  extra_physics=physics)
                if evolving_tide:
                    apply_stellar_tide(0)
                    backend.remote[0].tides_viscous_time_scale_prescription = 1
            if old is not None and backend is not None:
                backend.restore(cp["state"])
                wd_require_state(backend.checkpoint(), cp["state"], n0, full=True)
                if cp["anchor"] is not None:
                    prior = cp["anchor"]
                    backend.restore(prior["state"])
                    reproduced, event = advance_one(prior["time_yr"], actual)
                    if event or reproduced != actual:
                        raise RuntimeError("The checkpoint replay reached an unexpected event.")
                    wd_require_state(backend.checkpoint(), cp["state"], n0, full=True)
                fh_json(folder/"restart_check.json", {"status": "passed_full_state_restart",
                          "index": index, "one_interval_replayed": cp["anchor"] is not None})
            elif old is not None:
                # A terminal event can fire again if restored. Its flushed row is sufficient.
                fh_json(folder/"restart_check.json", {"status": "verified_terminal_checkpoint",
                          "index": index, "one_interval_replayed": False})
            for index in range(begin, len(targets)):
                event = ""
                save_due = (index == 0 or index == len(targets)-1
                            or perf_counter()-last_checkpoint >= limits["checkpoint_seconds"])
                anchor = ({"time_yr": actual, "state": backend.checkpoint()}
                          if save_due and index else None)
                if index:
                    actual, event = advance_one(actual, targets[index])
                mdot, _ = interval_rates(interval_index(actual))
                row = fh_row(backend.snapshot(), actual, driver, cfg, physics, mdot, limits)
                mass_error = abs(row["mass_MSun"]/row["driver_mass_MSun"]-1)
                radius_error = abs(row["radius_AU"]-row["driver_radius_AU"])
                if mass_error > limits["mass_transfer_relative"] or radius_error > (
                        limits["radius_transfer_absolute_AU"]
                        + limits["radius_transfer_relative"]*row["driver_radius_AU"]):
                    raise RuntimeError("Stellar mass/radius transfer exceeded the working limits.")
                if table is None:
                    table = np.lib.format.open_memmap(folder/"history.npy", mode="w+",
                        dtype=[(key, "f8") for key in row], shape=(capacity,))
                size = index+1
                phase_change = not phase_records or row["stellar_type"] != phase_records[-1]["stellar_type"]
                if phase_change:
                    phase_records.append({"time_yr": actual, "age_Myr": row["age_Myr"],
                                          "stellar_type": row["stellar_type"], "index": index})
                if limits["roche_enabled"] and roche_record is None and row["roche_clearance_AU"] <= 0:
                    roche_record = {"time_yr": actual, "age_Myr": row["age_Myr"], "index": index,
                                    "label": "first_sampled_roche_crossing"}
                if event:
                    outcome = event.replace("instability_criterion", "dynamical_instability")
                    event_origin = "solver_root"
                else:
                    sampled = fh_sampled_stop(row, cfg, physics, limits)
                    if sampled:
                        outcome, event_origin = sampled, "sampled_guard"
                fh_update_extrema(extrema, row, physics["companion"])
                if (index % history_stride == 0 or phase_change
                        or index == len(targets) - 1 or outcome != "completed_window"):
                    table[saved_size] = tuple(row[key] for key in table.dtype.names)
                    saved_size += 1
                if save_due or outcome != "completed_window":
                    checkpoint()
                if perf_counter()-last_print >= limits["progress_seconds"]:
                    print(f"{folder.name}: {size:,}/{len(targets):,} states, "
                          f"{(wall_before+perf_counter()-started)/60:.1f} min", flush=True)
                    last_print = perf_counter()
                if outcome != "completed_window":
                    break
            table.flush()
            summary = dict(run=folder.name, outcome=outcome, n_states=size, n_saved=saved_size,
                           history_stride=history_stride, end_time_yr=row["time_yr"],
                           end_age_Myr=row["age_Myr"], final_a_in_AU=row["a_in_AU"],
                           final_e_in=row["e_in"], final_a_out_AU=row["a_out_AU"] if physics["companion"] else None,
                           final_e_out=row["e_out"] if physics["companion"] else None, **extrema)
            if not physics["companion"]:
                summary.update(max_adiabatic_ratio_out=None, max_sampled_double_averaging_ratio=None,
                               minimum_sampled_hierarchy_ratio=None)
            wd_indices = np.flatnonzero(np.isin(driver["stellar_type"], (10, 11, 12)))
            first_wd = int(wd_indices[0]) if len(wd_indices) else None
            domain_exits = ("dynamical_instability", "instability_criterion", "orbit_crossing", "secular_breakdown",
                            "sampled_adiabatic_model_limit", "sampled_spin_model_limit",
                            "sampled_double_averaging_limit")
            summary.update(source_first_WD_age_Myr=(float(driver["age_Myr"][first_wd])
                                                      if first_wd is not None else None),
                source_wd_bracket_start_age_Myr=(float(driver["age_Myr"][first_wd-1])
                                                 if first_wd is not None and first_wd > 0 else None),
                within_model_domain=not any(k in outcome for k in domain_exits),
                event_origin=event_origin, requested_end_age_Myr=float(driver["age_Myr"][-1]),
                pilot=pilot,
                reached_wd_endpoint=bool(not pilot and first_wd is not None
                                         and outcome == "completed_window" and actual == targets[-1]),
                reached_requested_end=bool(outcome == "completed_window" and actual == targets[-1]),
                phase_records=phase_records, sampled_roche_crossing=roche_record,
                wall_seconds=wall_before+perf_counter()-started)
            fh_json(folder/"summary.json", summary)
            fh_json(folder/"run.json", {"status": "completed", "inputs": inputs,
                    "n_states": size, "n_saved": saved_size, "outcome": outcome,
                    "wall_seconds": summary["wall_seconds"]})
    except BaseException as exc:
        # Keep the last coherent checkpoint. Rows after it are overwritten only on explicit resume.
        if table is not None:
            table.flush()
        fh_json(folder/"run.json", {"status": "interrupted", "inputs": inputs,
                "n_states": size, "error": repr(exc),
                "wall_seconds": wall_before+perf_counter()-started})
        raise
    finally:
        if backend is not None:
            backend.stop()
    return fh_read_case(folder)


# Read SeBa histories and construct nested integration grids.
FH_RSUN_AU = 0.004650467260962158
FH_WD_TYPES = (10, 11, 12)


def fh_validate_track(track):
    keys = ("age_Myr", "mass_MSun", "radius_AU", "stellar_type")
    age = np.asarray(track["age_Myr"])
    if (age.ndim != 1 or len(age) < 2 or np.any(np.diff(age) <= 0)
        or any(np.shape(track[k]) != age.shape or not np.isfinite(track[k]).all() for k in keys)
        or np.any(track["mass_MSun"] <= 0) or np.any(track["radius_AU"] <= 0)
        or np.any(track["stellar_type"] != np.floor(track["stellar_type"]))):
        raise ValueError("Invalid stellar ages, masses, radii or phase labels.")
    if "time_yr" in track:
        expected = (age-age[0])*1e6
        if not np.allclose(track["time_yr"], expected, rtol=0, atol=1e-6):
            raise ValueError("The local clock does not match the stellar ages.")
    return track


def fh_sample_track(track, ages_Myr):
    # Radius and mass are linear between saved states; a stellar type is a label.
    age = np.asarray(ages_Myr, dtype=float)
    source = track["age_Myr"]
    if not np.isfinite(age).all() or np.any(age < source[0]) or np.any(age > source[-1]):
        raise ValueError("The requested ages are outside the saved stellar coverage.")
    i = np.searchsorted(source, age, side="right")-1
    return dict(mass_MSun=np.interp(age, source, track["mass_MSun"]),
                radius_AU=np.interp(age, source, track["radius_AU"]),
                stellar_type=track["stellar_type"][i],
                **{key: track[key][i] for key in STELLAR_STRUCTURE_FIELDS if key in track})


def fh_crop_track(track, start_age_Myr=None, end_age_Myr=None):
    age = track["age_Myr"]
    start = float(age[0] if start_age_Myr is None else start_age_Myr)
    end = float(age[-1] if end_age_Myr is None else end_age_Myr)
    if not age[0] <= start < end <= age[-1]:
        raise ValueError("The common window must lie inside the saved source; no extrapolation.")
    new_age = np.unique(np.r_[start, age[(age > start) & (age < end)], end])
    result = {k: v for k, v in track.items() if not isinstance(v, np.ndarray)}
    result.update(fh_sample_track(track, new_age))
    result.update(age_Myr=new_age, time_yr=(new_age-start)*1e6)
    return fh_validate_track(result)


def fh_coupling_counts(track, max_step_yr, max_mass_fraction=0.002,
                       max_radius_fraction=0.01, subdivisions=1):
    fh_validate_track(track)
    if (not np.isfinite([max_step_yr, max_mass_fraction, max_radius_fraction]).all()
        or min(max_step_yr, max_mass_fraction, max_radius_fraction) <= 0
        or int(subdivisions) != subdivisions or subdivisions < 1):
        raise ValueError("Use positive step limits and an integer subdivision factor.")
    t = track["time_yr"]
    counts = np.maximum(1., np.ceil(np.diff(t)/max_step_yr))
    for key, fraction in (("mass_MSun", max_mass_fraction), ("radius_AU", max_radius_fraction)):
        v = track[key]
        counts = np.maximum(counts, np.ceil(abs(np.diff(v))/np.minimum(v[:-1], v[1:])/fraction))
    return counts, int(subdivisions)


def fh_coupling_grid(track, max_step_yr, max_mass_fraction=0.002,
                     max_radius_fraction=0.01, subdivisions=1, max_outputs=5_000_000):
    counts, subdivisions = fh_coupling_counts(
        track, max_step_yr, max_mass_fraction, max_radius_fraction, subdivisions)
    count = int(counts.sum()) * subdivisions + 1
    if count > max_outputs:
        raise ValueError(f"This schedule needs {count:,} outputs, above max_outputs={max_outputs:,}. Review its physical step caps and runtime before raising the budget.")
    counts = counts.astype(np.int64)*int(subdivisions)
    t = track["time_yr"]
    grid = np.empty(count, dtype=float)
    grid[0] = t[0]
    # Reproduce linspace(left, right, n+1)[1:] in blocks. Most of the
    # 1.7 million driver intervals need only one or two points.
    offsets = np.r_[0, np.cumsum(counts)]
    for start in range(1, count, 65_536):
        stop = min(start + 65_536, count)
        position = np.arange(start, stop)
        interval = np.searchsorted(offsets, position, side="left") - 1
        local = position - offsets[interval]
        n = counts[interval]
        left, right = t[interval], t[interval + 1]
        step = (right - left) / n
        grid[start:stop] = np.where(local == n, right, left + local * step)
    if np.any(np.diff(grid) <= 0):
        raise ValueError("The requested schedule is finer than the available clock precision.")
    return grid
