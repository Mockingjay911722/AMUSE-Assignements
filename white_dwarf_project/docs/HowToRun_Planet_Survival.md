# Reproducing the white-dwarf planet-survival campaign

Start here to regenerate the data for the 72-case survey and its 56 follow-up configurations. SeBa evolves one star. Each SecularMultiple job follows a planet, with or without a brown dwarf, through that shared stellar history. The reference configuration schedules **214 independent integrations**. You can also run one case, change the model, or choose how much of each orbit to save.

## First run

Open a terminal in **`white_dwarf_project/`**. If you have an AMUSE environment, activate it. Installation instructions for a new machine are below.

```bash
python scripts/survival.py plan
python scripts/survival.py pilot --suite survey --case a1_isolated a1_i60 --level n1 --age-myr 1 --jobs 2
```

The first command lists the seven suites. The second generates the first 1 Myr of stellar evolution from scratch and runs two orbits concurrently. With the reference settings, each finishes with **78 integrated states and 11 saved rows**. The pilot goes in `data/survival_campaign/pilots/1.0Myr/`. A pilot that reaches 1 Myr has completed its test window.

For a small production example, the 1 au, 90° case reaches contact early. The following sequence generates the full reference star once, runs both resolutions of this case and compares them to the original CSV.

```bash
python scripts/survival.py prepare
python scripts/survival.py run --suite survey --case a1_i90 --jobs 2
python scripts/survival.py collect --compare
```

Expect two completed jobs and 212 pending. Both historical runs reached contact near **0.2329511 Myr**. The full campaign uses the same commands without a case filter.

```bash
python scripts/survival.py run --jobs auto
python scripts/survival.py collect --compare
```

`run` prepares the star automatically if needed. By default, every command uses `data/survival_campaign/` relative to the project. To use another disk, pass the same `--output /path/to/campaign` after each command. Use `python scripts/survival.py run --help` for the available run options.

## Installing on another machine

Use **Python 3.11 or later**, NumPy, the AMUSE framework, and working **SeBa** and **SecularMultiple** community workers. The simulation uses the `sockets` channel by default. An AMUSE environment already running notebooks 01 and 06 should have both workers. The basic plan and collection commands use only the Python standard library. A storage estimate with `plan --driver` also uses NumPy.

For a new installation on Linux, macOS, or Windows under WSL, start in a directory separate from the course project. The AMUSE installer detects your operating system and identifies compiler, MPI and library prerequisites. One way to select the release used for the reference work is below. Use the executable name for your installed Python 3.11+ version.

```bash
git clone --depth 1 --branch v2026.3.0 https://github.com/amusecode/amuse.git
cd amuse
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip wheel
./setup
```

Follow the dependency instructions printed by `./setup`. Once the framework and the two codes are enabled for installation, install them in the active environment.

```bash
./setup install amuse-framework
./setup install amuse-seba amuse-secularmultiple
```

Return to `white_dwarf_project/` and run the short pilot above. This exercises both compiled workers. The original orbital run recorded AMUSE 2026.3.0 and NumPy 2.5.3. The native checks for this package used Python 3.12.14, AMUSE 2026.3.0 and NumPy 2.3.5 on Linux. For Conda, platform-specific setup, or a different MPI installation, use the upstream instructions. [AMUSE installation guide](https://amuse.readthedocs.io/en/latest/install/installing.html)

## Which runs are included?

`configs/survival_reference.toml` defines all initial conditions and the numerical levels scheduled for each case. The six follow-up suites explore alternatives inspired by KMT-2020-BLG-0414. The adaptive suite replays the branches actually run, including all three tight refinements. Three candidates left pending in the original campaign are excluded.

| Suite | Configurations | n1 | n2 | n2_tight | Save every |
| --- | ---: | ---: | ---: | ---: | ---: |
| `survey` · planet distance and mutual inclination | 72 | 72 | 72 | 0 | 8 states |
| `kmt` · companion distance and tidal controls | 20 | 20 | 0 | 0 | 8 states |
| `kmt-finite` · 30 and 100 yr stellar viscous timescales | 4 | 4 | 0 | 0 | 8 states |
| `kmt-adaptive` · selected weak-tide branches | 4 | 4 | 3 | 3 | 8 states |
| `kmt-inclinations` · intermediate inclinations and isolated control | 4 | 4 | 3 | 2 | 8 states |
| `kmt-eccentric` · eccentric companion and periapse angle | 16 | 16 | 2 | 0 | 16 states |
| `kmt-envelope` · evolving-envelope stellar tides | 8 | 8 | 1 | 0 | 16 states |
| **Total** | **128** | **128** | **81** | **5** | |

Use `plan --list` to see exact case labels and scheduled levels. Filters can select several names.

```bash
python scripts/survival.py plan --suite kmt-eccentric --level n2 --list
python scripts/survival.py run --suite kmt-eccentric --case a1_i40_outer15_e0.3_w90_star_tv300 --level n2 --jobs 1
```

### Three time intervals you can change

| What you want to change | Setting | Meaning |
| --- | --- | --- |
| How often the star is sampled | `[stellar]` in TOML | The main sequence uses the phase-specific base spacing divided by 128. After the first saved AGB state, the reference refinement is 2048. SeBa continues through the first 100 WD years at 0.5 yr intervals. |
| How often the star and orbit are coupled | `n1`, `n2` or custom levels, plus `[solver]` | Every stellar knot is included. Additional subdivisions limit elapsed time and fractional stellar mass/radius changes. `n2` halves each `n1` interval. A `_tight` suffix also divides both solver tolerances by 10 by default. |
| How often an orbital state is written to disk | `--history-stride N` or `history_stride` per suite | Write row zero, every Nth coupling row, the first sample of each stellar phase, and the terminal row. N = 1 gives a dense history. Adaptive coupling means the elapsed years between saved rows can vary. |

For example, this saves a dense 1 au, 60° history using the already prepared star. Use a separate output directory when changing the save cadence for an existing job.

```bash
python scripts/survival.py run --suite survey --case a1_i60 --level n2 --history-stride 1 \
  --driver data/survival_campaign/stellar_driver.npz --output data/survival_dense --jobs 1
```

The solver advances, checks stopping conditions and updates extrema at every coupling step for every stride. The smallest clearance or largest eccentricity can occur between saved rows. Use `summary.json` or the collected CSV for those measurements. Checkpoint frequency is elapsed **wall time** and defaults to 60 seconds.

### Changing the physics or case grid

Copy `configs/survival_reference.toml`, edit the copy and pass `--config` after the subcommand. `[stellar]` sets the initial stellar mass, metallicity, sampling and endpoint. `[orbit]` supplies the common planet and companion masses, radii, orbital elements, spin rates and tidal parameters. `[physics]` sets relativistic precession, octupole gravity and the stellar wind spin term. The reference mass is 1.2 solar masses at Z = 0.02, with a 1.9 Earth-mass planet and a 27 Jupiter-mass brown dwarf. The star and planet start with the same physical spin rates across the main survey's distance grid.

The groups in each suite take the Cartesian product of their arrays. Add `isolated = true` to include a companion-off control. `levels` applies to the entire suite and `refine` adds specific case/level pairs. For example, append this to the **copy** to test a 1.7 au planet at a new inclination:

```toml
[suites.my_grid]
history_stride = 4
levels = ["n1", "n2"]

[[suites.my_grid.groups]]
a_in_AU = [1.7]
inclinations_deg = [35]
isolated = true
```

```bash
python scripts/survival.py plan --config configs/my_campaign.toml --suite my_grid --list
python scripts/survival.py run --config configs/my_campaign.toml --suite my_grid --output data/my_grid \
  --driver data/survival_campaign/stellar_driver.npz --jobs 2
python scripts/survival.py collect --output data/my_grid
```

Reuse the prepared driver when the stellar settings match. If you change the star, omit `--driver` to generate the new star in the new output directory. The orbital reference tolerances are `rtol = 1e-14` and an absolute eccentricity-vector tolerance of `1e-16`. The reference `n1` grid caps fractional stellar mass and radius changes at 0.002 and 0.01, respectively. The maximum interval also depends on one quarter of the initial Lidov–Kozai timescale. The TOML records explicit settings for every group.

## Time, parallelism and disk space

| Work | Wall time or cost | Basis |
| --- | ---: | --- |
| Generate the first 500 Myr of the reference star | 16 s | Measured locally with native SeBa. All 12,801 rows agreed exactly with the machine export in all five shared columns. |
| Prepare the full reference star once | Budget roughly an hour for an initial attempt | The original AGB-2048 stellar run took about 49 min. Scaling the local short measurement to 1,668,280 samples gives about 35 min. The new preparation also reads the envelope. Hardware and stellar phase affect this estimate. |
| Run 1 Myr of the isolated and 60° orbits, `n1` | About 2 s for the pair | Measured locally with two native orbit workers, including their startup. Stellar preparation and software installation are separate. |
| Original notebook 06, 1 au, 60°, dense fine run | **4 h 21 min 35 s** | Measured wall time from the uploaded `run.json`, 1,828,043 integrated states before contact. |
| Campaign 1 au, 60° run | About **2 h 11 min** at `n1`, **4 h 22 min** at `n2` | Estimated from 914,058 and 1,828,041 archived states using the original rate of 116.47 states/s. |
| Campaign 2 au, isolated, `n2` | About **9 h 20 min** | Same estimate for 3,914,275 states through the WD endpoint. |
| All 144 main-survey jobs | About **698 serial job-hours** | Same estimate for 292,823,908 archived states. Eight equally fast, continuously occupied workers would take about 87 hours. |

Many highly inclined cases reach an early event. Some wider triples encounter a model limit at the initial state. Neither gives a representative throughput for a full run. To benchmark your machine, time **the same pair** in fresh output directories at two worker counts. The command-line batch time includes worker startup, integration and file I/O after the shared star is ready.

```bash
python scripts/survival.py pilot --suite survey --case a1_isolated a1_i60 --level n1 --age-myr 50 --jobs 1 --output data/timing_1
python scripts/survival.py pilot --suite survey --case a1_isolated a1_i60 --level n1 --age-myr 50 --jobs 2 --output data/timing_2
```

Each concurrent job has its own Python process and AMUSE worker. Jobs across numerical levels are independent. `--jobs auto` uses the detected CPU allocation and available memory on Linux, budgeting 2 GiB per job plus 1 GiB for the host. When available memory cannot be determined, auto starts at most two jobs. `--jobs N` overrides the estimate. Choose N to fit your machine or scheduler allocation. A compute node, a job scheduler and an SSH connection are all optional. A single job uses one orbital worker.

With the reference driver and original strides, the full 214-job plan allows up to **35.54 GiB** of trajectory arrays. Allow about **40 GiB** for these histories, the 115 MiB prepared star, metadata and logs. The 208 archived numerical records actually report 439,996,012 integrated states and 51,021,828 saved rows, about 17.87 GiB of valid history data. Each history file is sized for its planned capacity at startup, and filesystem allocation for unused rows varies. Dense output increases this budget substantially. Get the estimate for your selection first:

```bash
python scripts/survival.py plan --driver data/survival_campaign/stellar_driver.npz
python scripts/survival.py plan --suite survey --case a1_i60 --level n2 --history-stride 1 \
  --driver data/survival_campaign/stellar_driver.npz --list
```

## Results, restart and comparison

| Output relative to `--output` | What it contains |
| --- | --- |
| `stellar_driver.npz` | SeBa ages, masses, radii, stellar types, requested times, luminosities, convective-envelope masses/radii, and preparation metadata. Created when generating the star locally. |
| `campaign.json` | The resolved TOML configuration and SHA-256 of the input driver. |
| `runs/<suite>/<case>/<level>/` | `job.json` for the case definition, `worker.log` for progress, `checkpoint.json` for restart, `run.json` for status and timing, `history.npy` for the saved states, and `summary.json` on completion. |
| `planet_survival_results.csv` | Created by `collect`. One row per completed production job, including cross-level outcomes and contact-age shifts where available. |
| `reference_comparison.csv` | Created by `collect --compare`. Per-job status, outcome agreement and signed numerical differences (rerun minus archive). |

Repeat the same run command after interruption. Successful jobs are reused. A nonterminal orbital checkpoint restores the orbital vectors, stellar state and spins, and replays its preceding interval when an anchor is available. A failed startup before the first checkpoint is retained under `.attempts/` and retried. Stellar preparation publishes its NPZ after the entire pass and starts again from age zero if interrupted. Use a new output directory for a changed campaign. You can change `--jobs` or select additional cases from the same campaign when you resume.

Read just the **valid prefix** of a finished history. The memory map gives access to named columns without loading the entire file.

```python
import json
from pathlib import Path
import numpy as np

case = Path("data/survival_campaign/runs/survey/a1_i90/n2")
run = json.loads((case / "run.json").read_text())
history = np.load(case / "history.npy", mmap_mode="r", allow_pickle=False)[:run["n_saved"]]
print(history.dtype.names)
print(history[["age_Myr", "a_in_AU", "e_in", "clearance_AU"]][-1])
```

A surface-contact result means the planetary pericentre reached the sum of the stellar and planetary radii. A WD endpoint means the orbit reached the end of the supplied stellar history. A model stop, such as the outer period exceeding 1% of the estimated Lidov–Kozai timescale, leaves subsequent survival undetermined. In the reference data the first saved WD occurs at **6514.110102252153 Myr** and the endpoint is **6514.110202252154 Myr**.

### Using the original files

The package includes the original `reference/planet_survival_results.csv` and notebook `reference/run.json` for comparison. The campaign generates its scientific inputs from scratch. If you have the optional machine export, `driver.npz`, it gives two additional shortcuts:

```bash
python scripts/survival.py prepare --verify-against /path/to/driver.npz
python scripts/survival.py run --suite survey --case a1_i90 --driver /path/to/driver.npz \
  --output data/archive_check --jobs 2
python scripts/survival.py collect --output data/archive_check --compare
```

The first command checks all five common stellar columns across the whole track. A generated NPZ records its stellar settings. A separate five-column file of unknown origin can establish mass and age coverage, but its metallicity and sampling provenance must be supplied by its author. The time tolerances are 0.001 yr and 1e-9 Myr, stellar types must match exactly, and mass/radius use 1e-9 relative tolerances with absolute floors of 1e-12 solar masses and 1e-14 au. The second command can run immediately with the original five-column file. Hurley tides require the envelope fields obtained by `prepare`. For a campaign that includes that suite, use the enriched driver from the outset. The original five-column file has SHA-256 `340e5e0c992df4929dbc097d83a5a297e9ddb00899d4ff6eaeef096592f1049b`.

The dense notebook 06 calculation used a common time grid constructed from two stellar refinements. The survey used the single reference stellar track and its case-specific coupling grid. For the 1 au, 60° fine run, those grids contain 3,914,399 and 3,914,275 planned states, respectively. The campaign comparison target is the historical survey CSV. It contains 208 numerical records and six rows preserving only the outcome. A fresh completed run supplies numerical measurements for each of those six, and the comparison leaves the unavailable historical differences blank.

## Verification

```bash
python -m unittest discover -s tests -v
SURVIVAL_REFERENCE_DRIVER=/path/to/driver.npz python -m unittest discover -s tests -v
```

The tests cover the exact set of archived job identifiers, reference physics, historical save cadences and the CLI. Supplying the optional machine driver adds a full-track coupling-count check. With the uploaded full stellar driver, 30 short production runs reproduced the archived outcome, end age, minimum clearance, peak eccentricity and row counts exactly. They covered 14 early contacts across seven planet distances and 16 initial double-averaging model stops at 5 au, including both coupling levels. At `n2`, strides 1, 8, 16 and 1,000 saved 37, 6, 4 and 2 rows from the same 37 integrated states. Native tests also exercised the Hurley model, checkpoint recovery and the first 500 Myr of the reference SeBa evolution. A separate coarse-grid example generated a star from age zero through its WD tail, evolved an isolated planet to that endpoint and produced a completed results CSV. Its grid and outer reference distance were changed to make this an end-to-end software check. The complete 214-job campaign and full-resolution SeBa regeneration were not rerun for this release.
