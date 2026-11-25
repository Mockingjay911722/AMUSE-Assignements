"""Regression checks for the published configuration. No AMUSE worker is needed."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from campaign import build_jobs, load_config, select_jobs, settings_for
from orbit_worker import load_driver
from survey_orbits import fh_coupling_counts


def identity(row):
    return row["suite"], row["case"], row["data_level"]


class ReferenceCampaign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = load_config()
        cls.jobs = build_jobs(cls.config)
        cls.by_key = {identity(job): job for job in cls.jobs}
        archive_path = ROOT / "reference/planet_survival_results.csv"
        with archive_path.open(newline="", encoding="utf-8") as f:
            cls.archive = list(csv.DictReader(f))

    def test_exact_reference_schedule(self):
        self.assertCountEqual(map(identity, self.jobs), map(identity, self.archive))
        self.assertEqual(len(self.by_key), 214)

    def test_historical_save_cadence(self):
        for row in self.archive:
            with self.subTest(job=identity(row)):
                stride = self.by_key[identity(row)]["history_stride"]
                wide_stride = row["suite"] in ("kmt-eccentric", "kmt-envelope")
                self.assertEqual(stride, 16 if wide_stride else 8)
                if row["integrated_states"] and row["saved_history_rows"]:
                    # Five stellar type transitions and one terminal row may be off-stride.
                    regular = 1 + (int(row["integrated_states"]) - 1) // stride
                    self.assertIn(int(row["saved_history_rows"]) - regular, range(7))

    def test_save_override_preserves_physics_and_coupling(self):
        changed = select_jobs(self.jobs, history_stride=3)
        self.assertEqual(changed, [{**job, "history_stride": 3} for job in self.jobs])
        self.assertEqual({job["history_stride"] for job in self.jobs}, {8, 16})
        with self.assertRaises(ValueError):
            select_jobs(self.jobs, history_stride=0)

    def test_cli_plan_outside_project(self):
        script = ROOT / "scripts/survival.py"
        command = [sys.executable, str(script), "plan", "--suite", "survey", "--case", "a1_i60",
                   "--level", "n2", "--history-stride", "3", "--list"]
        with tempfile.TemporaryDirectory() as cwd:
            result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertRegex(result.stdout, r"(?m)^survey\s+a1_i60\s+n2\s+stride=3\b")
            self.assertEqual(list(Path(cwd).iterdir()), [])

    def test_shared_model_matches_original_run(self):
        record = json.loads((ROOT / "reference/run.json").read_text(encoding="utf-8"))
        settings = settings_for(self.config)
        for name in ("cfg", "physics", "rtol", "atol", "worker_options", "limits"):
            with self.subTest(setting=name):
                self.assertEqual(settings[name], record["inputs"][name])
        self.assertEqual(settings["max_mass_fraction"], 0.002)
        self.assertEqual(settings["max_radius_fraction"], 0.01)
        self.assertEqual(settings["tight_factor"], 10)

    def test_distinct_tide_models_and_orientations(self):
        off = self.by_key["kmt", "a1_i0_outer15_friction_off", "n1"]["spec"]["cfg"]
        self.assertFalse(off["tides"])
        self.assertTrue(off["tidal_bulges"] and off["rotation_bulges"])
        key = ("kmt-adaptive", "a1_i0_outer15_planet_tides_only", "n2")
        planet = self.by_key[key]["spec"]["cfg"]
        self.assertTrue(planet["tides"])
        self.assertFalse(planet["tides_star"])
        isolated = self.by_key["kmt-envelope", "a1_isolated_hurley", "n2"]["spec"]
        self.assertIsNone(isolated["inclination_deg"])
        self.assertEqual(isolated["cfg"]["stellar_tide_prescription"], 1)
        key = ("kmt-eccentric", "a1_i40_outer15_e0.3_w90_star_tv300", "n2")
        eccentric = self.by_key[key]["spec"]["cfg"]
        outer = (eccentric["a_out_AU"], eccentric["e_out"], eccentric["omega_out_deg"])
        self.assertEqual(outer, (15, 0.3, 90))

    @unittest.skipUnless(os.environ.get("SURVIVAL_REFERENCE_DRIVER"),
                         "set SURVIVAL_REFERENCE_DRIVER to the original driver.npz")
    def test_original_driver_target_count(self):
        driver = load_driver(os.environ["SURVIVAL_REFERENCE_DRIVER"])
        job = self.by_key["survey", "a1_i60", "n2"]
        settings = settings_for(self.config)
        count, level = fh_coupling_counts(driver, job["spec"]["step_yr"],
                                          max_mass_fraction=settings["max_mass_fraction"],
                                          max_radius_fraction=settings["max_radius_fraction"],
                                          subdivisions=job["level"])
        self.assertEqual(len(driver["time_yr"]), 1_668_280)
        self.assertEqual(1 + int(count.sum()) * level, 3_914_275)


if __name__ == "__main__":
    unittest.main()
