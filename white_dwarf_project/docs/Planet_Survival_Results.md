## Full-history planet survival results

The 72-case survey varies the planet's starting distance and the angle between the two orbits (planet and brown dwarf). The 56 follow-up cases use smaller grids to see what changes when we vary the tides or the brown dwarf's orbit. These are full-history runs (like in notebook 06) with the same SeBa stellar track for a 1.2 M☉ star at Z = 0.02. The planet is about 1.9 Earth masses and the companion, where included, is about 27 Jupiter masses. The `n2` runs update the orbit with the changing stellar input twice as often as `n1`. Some interesting results are:

- **The tide model matters even in isolation.** At 1 au, the planet reaches the white dwarf with a constant stellar viscous timescale of 300 yr and reaches the stellar surface with evolving-envelope Hurley tides (These, as I understand them, are differential tides, like the bulging and receding that happens to the earth, most noticably the ocean in the equator due to the influence of the moon). Both outcomes agree at the two coupling grids.
- **The brown dwarf can change the outcome at fixed tides.** For a nearly circular companion starting at 15 au and the same 300 yr setting, the 40° case reaches the endpoint and the 50° case reaches contact. Both outcomes agree across the tested coupling grids and solver settings.
- **Some outcomes are sensitive to the time step.** At 40° and outer eccentricity 0.3, the case with initial outer ω = 0° reaches contact at `n1` and the WD endpoint at `n2`. The ω = 90° case reaches the endpoint on both grids.

**The tables below.** WD means the run reached stellar age 6514.110202 Myr, about 100 years after the first saved white dwarf state. Contact means the planet's pericentre, `q = a(1 − e)`, reached the combined stellar and planetary radii. The minimum sampled clearance is the smallest value of `q − (Rstar + Rplanet)` checked at a coupling step. “Disagree” flags a change in outcome across the available settings. “—” marks a not existent entry, mostly appears due to incompatible settings.

**Tight runs**  The tight runs use the same `n2` coupling grid and divide two SecularMultiple solver tolerances by 10. The notebook 06 reference settings give the following values.

| Solver setting | Standard (`n1` and `n2`) | Tight |
| --- | --- | --- |
| Relative tolerance (`rtol`) | `1e-14` | `1e-15` |
| Absolute tolerance on eccentricity-vector components (`atol`) | `1e-16` | `1e-17` |

The relative tolerance scales the solver's error target with the size of the quantity it is integrating. The absolute tolerance sets the error scale for eccentricity-vector components near zero.

**What does “Model stop” in the tables mean?** SecularMultiple averages the gravitational interaction over both orbits and follows the slower changes in their shapes and orientations. This works best when the brown dwarf completes many orbits in one Lidov–Kozai timescale, `t_LK`, which estimates how quickly the companion can change the planet's orbit. At each coupling sample, our runner compares the brown dwarf's orbital period, `P_out`, with `t_LK`. It stops at the first sample where `P_out / t_LK > 0.01`, equivalent to fewer than 100 outer orbits per estimated Lidov–Kozai timescale. This is our cutoff for trusting the averaged model. We probably should consider these runs undetermined.

### Main survey

The companions start at 20 au with eccentricity 0.001. These runs use the reference tides with a stellar viscous timescale of 1 yr. All 72 cases have both coupling levels.

| Initial planet orbit [au] | Isolated | 0° | 30° | 45° | 60° | 70° | 80° | 90° | 120° |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.8 | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact |
| 1 | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact |
| 1.2 | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact |
| 1.5 | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact | Contact |
| 2 | WD | WD | WD | **Disagree** | Contact | Contact | Contact | Contact | Contact |
| 2.5 | WD | Model stop | Model stop | Model stop | Contact | Contact | Contact | Contact | Contact |
| 3.5 | WD | Model stop | Model stop | Model stop | **Disagree** | Contact | Contact | Contact | Model stop |
| 5 | WD | Model stop | Model stop | Model stop | Model stop | Model stop | Model stop | Model stop | Model stop |

There are 49 contacts, 6 cases reaching the WD endpoint, 15 model stops and 2 disagreements. The 2 au, 45° case goes from contact to WD, and the 3.5 au, 60° case goes from contact to a model stop as the coupling steps get smaller.

<details>
<summary>Numerical results for all 72 main-survey cases</summary>

These values come from `n2`. Contact age shift is the absolute difference between the two coupling levels where both reached contact.

| Initial planet orbit [au] | Configuration | `n1` | `n2` | Peak inner eccentricity | Min. clearance [au] | Stop age [Myr] | Contact age shift [yr] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.8 | Isolated | Contact | Contact | 0.001000 | ≈0 | 6375.051664 | 0.013 |
| 0.8 | 0° | Contact | Contact | 0.001000 | ≈0 | 6375.051665 | 0.005 |
| 0.8 | 30° | Contact | Contact | 0.001591 | ≈0 | 6375.051651 | 0.305 |
| 0.8 | 45° | Contact | Contact | 0.391699 | ≈0 | 6374.282536 | 18.658 |
| 0.8 | 60° | Contact | Contact | 0.756360 | ≈0 | 6366.099793 | 28162.766 |
| 0.8 | 70° | Contact | Contact | 0.892289 | ≈0 | 6340.035717 | 31536.137 |
| 0.8 | 80° | Contact | Contact | 0.973090 | ≈0 | 6023.618531 | 24032474.458 |
| 0.8 | 90° | Contact | Contact | 0.993293 | ≈0 | 0.324183 | 0.000 |
| 0.8 | 120° | Contact | Contact | 0.756368 | ≈0 | 6366.072947 | 32413.195 |
| 1 | Isolated | Contact | Contact | 0.001000 | ≈0 | 6376.352495 | 0.010 |
| 1 | 0° | Contact | Contact | 0.001000 | ≈0 | 6376.352496 | 0.003 |
| 1 | 30° | Contact | Contact | 0.001604 | ≈0 | 6376.352487 | 0.201 |
| 1 | 45° | Contact | Contact | 0.401632 | ≈0 | 6375.740135 | 1.309 |
| 1 | 60° | Contact | Contact | 0.760796 | ≈0 | 6369.836990 | 57393.220 |
| 1 | 70° | Contact | Contact | 0.895258 | ≈0 | 6351.617915 | 31372.706 |
| 1 | 80° | Contact | Contact | 0.973460 | ≈0 | 6152.141581 | 6938708.954 |
| 1 | 90° | Contact | Contact | 0.994639 | ≈0 | 0.232951 | 0.000 |
| 1 | 120° | Contact | Contact | 0.760825 | ≈0 | 6369.953090 | 11540.972 |
| 1.2 | Isolated | Contact | Contact | 0.001000 | ≈0 | 6377.142668 | 0.007 |
| 1.2 | 0° | Contact | Contact | 0.001000 | ≈0 | 6377.142669 | 0.002 |
| 1.2 | 30° | Contact | Contact | 0.001606 | ≈0 | 6377.142662 | 0.208 |
| 1.2 | 45° | Contact | Contact | 0.405218 | ≈0 | 6376.721447 | 20369.422 |
| 1.2 | 60° | Contact | Contact | 0.762343 | ≈0 | 6372.261547 | 170488.563 |
| 1.2 | 70° | Contact | Contact | 0.896314 | ≈0 | 6358.689523 | 28285.362 |
| 1.2 | 80° | Contact | Contact | 0.974057 | ≈0 | 6221.052754 | 3612883.763 |
| 1.2 | 90° | Contact | Contact | 0.995534 | ≈0 | 0.177583 | 0.000 |
| 1.2 | 120° | Contact | Contact | 0.762385 | ≈0 | 6372.207001 | 29411.064 |
| 1.5 | Isolated | Contact | Contact | 0.001000 | ≈0 | 6513.415203 | 0.030 |
| 1.5 | 0° | Contact | Contact | 0.001000 | ≈0 | 6513.415203 | 0.005 |
| 1.5 | 30° | Contact | Contact | 0.001604 | ≈0 | 6513.415200 | 0.148 |
| 1.5 | 45° | Contact | Contact | 0.407048 | ≈0 | 6377.481171 | 105231.383 |
| 1.5 | 60° | Contact | Contact | 0.763215 | ≈0 | 6374.440772 | 11813.765 |
| 1.5 | 70° | Contact | Contact | 0.896899 | ≈0 | 6365.082283 | 163398.524 |
| 1.5 | 80° | Contact | Contact | 0.974395 | ≈0 | 6273.083889 | 171189.656 |
| 1.5 | 90° | Contact | Contact | 0.996429 | ≈0 | 0.127276 | 0.000 |
| 1.5 | 120° | Contact | Contact | 0.763257 | ≈0 | 6374.552939 | 64832.586 |
| 2 | Isolated | WD | WD | 0.001000 | 1.44940 | 6514.110202 | — |
| 2 | 0° | WD | WD | 0.001000 | 1.44944 | 6514.110202 | — |
| 2 | 30° | WD | WD | 0.001596 | 1.44850 | 6514.110202 | — |
| 2 | 45° | Contact | WD | 0.407847 | 0.42133 | 6514.110202 | — |
| 2 | 60° | Contact | Contact | 0.763639 | ≈0 | 6376.334559 | 40346.770 |
| 2 | 70° | Contact | Contact | 0.897199 | ≈0 | 6370.263903 | 235467.213 |
| 2 | 80° | Contact | Contact | 0.974575 | ≈0 | 6316.628598 | 87778.882 |
| 2 | 90° | Contact | Contact | 0.997322 | ≈0 | 0.082774 | 0.000 |
| 2 | 120° | Contact | Contact | 0.763691 | ≈0 | 6376.317485 | 10344.879 |
| 2.5 | Isolated | WD | WD | 0.001000 | 2.06351 | 6514.110202 | — |
| 2.5 | 0° | Model stop | Model stop | 0.001000 | 2.06426 | 6513.909057 | — |
| 2.5 | 30° | Model stop | Model stop | 0.001588 | 2.06367 | 6513.909057 | — |
| 2.5 | 45° | Model stop | Model stop | 0.408261 | 0.84341 | 6513.912639 | — |
| 2.5 | 60° | Contact | Contact | 0.763772 | ≈0 | 6377.263083 | 183103.075 |
| 2.5 | 70° | Contact | Contact | 0.897296 | ≈0 | 6372.809636 | 20035.629 |
| 2.5 | 80° | Contact | Contact | 0.974647 | ≈0 | 6337.438732 | 309074.663 |
| 2.5 | 90° | Contact | Contact | 0.997858 | ≈0 | 0.059263 | 0.000 |
| 2.5 | 120° | Contact | Contact | 0.763831 | ≈0 | 6377.254877 | 1822.493 |
| 3.5 | Isolated | WD | WD | 0.001000 | 3.21200 | 6514.110202 | — |
| 3.5 | 0° | Model stop | Model stop | 0.001000 | 3.25137 | 6376.757031 | — |
| 3.5 | 30° | Model stop | Model stop | 0.001570 | 3.24936 | 6376.757031 | — |
| 3.5 | 45° | Model stop | Model stop | 0.408420 | 1.69912 | 6376.757812 | — |
| 3.5 | 60° | Contact | Model stop | 0.763874 | 0.16118 | 6513.318108 | — |
| 3.5 | 70° | Contact | Contact | 0.897389 | ≈0 | 6375.167330 | 1413.581 |
| 3.5 | 80° | Contact | Contact | 0.974727 | ≈0 | 6352.542678 | 87216.288 |
| 3.5 | 90° | Contact | Contact | 0.998470 | ≈0 | 0.035790 | 0.000 |
| 3.5 | 120° | Model stop | Model stop | 0.763944 | 0.17158 | 6512.718274 | — |
| 5 | Isolated | WD | WD | 0.001000 | 4.87984 | 6514.110202 | — |
| 5 | 0° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 30° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 45° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 60° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 70° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 80° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 90° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |
| 5 | 120° | Model stop | Model stop | 0.001000 | 4.98965 | 0.000000 | — |

</details>

### Isolated planet and inclination checks

The planet starts at 1 au. The triples here have a nearly circular 15 au companion orbit. The isolated Hurley row also appears in the full Hurley table below.

| Stellar tides | Configuration | Initial outer e | Initial outer ω | `n1` | `n2` | Tight | Result | Values from | Peak inner eccentricity | Min. clearance [au] | Stop age [Myr] | Final planet a [au] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 300 yr | Isolated | — | — | WD | WD | — | WD | n2 | 0.001000 | 0.17959 | 6514.110202 | 2.08798 |
| 300 yr | 30° | 0.001 | — | WD | — | — | WD | n1 | 0.001604 | 0.17957 | 6514.110202 | 2.08798 |
| 300 yr | 40° | 0.001 | — | WD | WD | WD | WD | n2 | 0.140023 | 0.12209 | 6514.110202 | 2.08054 |
| 300 yr | 50° | 0.001 | — | Contact | Contact | Contact | Contact | n2 | 0.556000 | ≈0 | 6376.288830 | — |
| Hurley | Isolated | — | — | Contact | Contact | — | Contact | n2 | 0.001000 | ≈0 | 6377.588562 | — |

The isolated Hurley contact ages differ by 0.00387 yr between `n1` and `n2`, about 1.4 days. The isolated 300 yr control has a minimum clearance of 0.17959 au at `n2`. A larger stellar viscous timescale gives weaker stellar friction at the same system state in the constant-tide prescription.

### Eccentric brown-dwarf orbits

The planet and brown dwarf start at 1 and 15 au. The stellar viscous timescale is fixed at 300 yr. The rows vary the mutual inclination, initial companion eccentricity and argument of periapse ω, which sets the direction of closest approach in the companion's orbital plane.

| Configuration | Initial outer e | Initial outer ω | `n1` | `n2` | Tight | Result | Values from | Peak inner eccentricity | Min. clearance [au] | Stop age [Myr] | Final planet a [au] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0° | 0.1 | 0° | WD | — | — | WD | n1 | 0.015726 | 0.16541 | 6514.110202 | 2.08763 |
| 30° | 0.1 | 0° | WD | — | — | WD | n1 | 0.016178 | 0.17375 | 6514.110202 | 2.08770 |
| 40° | 0.1 | 0° | WD | — | — | WD | n1 | 0.163007 | 0.02290 | 6514.110202 | 2.05574 |
| 50° | 0.1 | 0° | Contact | — | — | Contact | n1 | 0.565811 | ≈0 | 6376.292070 | — |
| 0° | 0.1 | 90° | WD | — | — | WD | n1 | 0.016786 | 0.16258 | 6514.110202 | 2.08755 |
| 30° | 0.1 | 90° | WD | — | — | WD | n1 | 0.015795 | 0.16515 | 6514.110202 | 2.08773 |
| 40° | 0.1 | 90° | WD | — | — | WD | n1 | 0.164203 | 0.09339 | 6514.110202 | 2.07454 |
| 50° | 0.1 | 90° | Contact | — | — | Contact | n1 | 0.565793 | ≈0 | 6376.288021 | — |
| 0° | 0.3 | 0° | WD | — | — | WD | n1 | 0.053745 | 0.12692 | 6514.110202 | 2.08358 |
| 30° | 0.3 | 0° | WD | — | — | WD | n1 | 0.056510 | 0.12351 | 6514.110202 | 2.08438 |
| 40° | 0.3 | 0° | Contact | WD | — | **Disagree** | n2 | 0.186545 | 0.07952 | 6514.110202 | 2.06513 |
| 50° | 0.3 | 0° | Contact | — | — | Contact | n1 | 0.609715 | ≈0 | 6376.220005 | — |
| 0° | 0.3 | 90° | WD | — | — | WD | n1 | 0.054766 | 0.12399 | 6514.110202 | 2.08348 |
| 30° | 0.3 | 90° | WD | — | — | WD | n1 | 0.050378 | 0.14340 | 6514.110202 | 2.08547 |
| 40° | 0.3 | 90° | WD | WD | — | WD | n2 | 0.190213 | 0.03757 | 6514.110202 | 2.05627 |
| 50° | 0.3 | 90° | Contact | — | — | Contact | n1 | 0.610149 | ≈0 | 6375.840883 | — |

### Evolving-envelope Hurley tides

The stellar tide responds to the changing stellar structure. Every case tested here reaches contact at `n1`. The isolated case also reaches contact at `n2`.

| Configuration | Initial outer e | Initial outer ω | `n1` | `n2` | Tight | Result | Values from | Peak inner eccentricity | Min. clearance [au] | Stop age [Myr] | Final planet a [au] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Isolated | — | — | Contact | Contact | — | Contact | n2 | 0.001000 | ≈0 | 6377.588562 | — |
| 0° | 0.001 | 0° | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6377.588558 | — |
| 40° | 0.001 | 0° | Contact | — | — | Contact | n1 | 0.140065 | ≈0 | 6377.559027 | — |
| 50° | 0.001 | 0° | Contact | — | — | Contact | n1 | 0.556034 | ≈0 | 6376.219062 | — |
| 0° | 0.3 | 0° | Contact | — | — | Contact | n1 | 0.053733 | ≈0 | 6377.579036 | — |
| 30° | 0.3 | 0° | Contact | — | — | Contact | n1 | 0.056492 | ≈0 | 6377.581702 | — |
| 40° | 0.3 | 0° | Contact | — | — | Contact | n1 | 0.185986 | ≈0 | 6377.548943 | — |
| 50° | 0.3 | 0° | Contact | — | — | Contact | n1 | 0.610009 | ≈0 | 6375.821245 | — |

### Other tidal controls

<details>
<summary>28 additional cases with weaker stellar tides, friction switched off, or planetary tides only</summary>

“Friction off” switches off damping on the star and planet. “Planet tides only” keeps damping on the planet. The numbers in the tides column give the adopted stellar viscous timescale. Three adaptive candidates were left pending in the campaign record, one at 60° with 1000 yr and two at 0° and 60° with 3000 yr.

| Initial planet a [au] | Initial companion a [au] | Configuration | Tides | `n1` | `n2` | Tight | Result | Values from | Peak inner eccentricity | Min. clearance [au] | Stop age [Myr] | Final planet a [au] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | — | Isolated | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6377.503859 | — |
| 1 | — | Isolated | Friction off | WD | — | — | WD | n1 | 0.001000 | 0.23064 | 6514.110202 | 2.22129 |
| 1.2 | — | Isolated | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6513.504235 | — |
| 1.2 | — | Isolated | Friction off | WD | — | — | WD | n1 | 0.001000 | 0.48231 | 6514.110202 | 2.66555 |
| 1.5 | — | Isolated | 10 yr | WD | — | — | WD | n1 | 0.001000 | 0.80912 | 6514.110202 | 3.08700 |
| 1.5 | — | Isolated | Friction off | WD | — | — | WD | n1 | 0.001000 | 0.85974 | 6514.110202 | 3.33193 |
| 1 | 15 | 0° | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6377.503861 | — |
| 1 | 15 | 60° | 10 yr | Contact | — | — | Contact | n1 | 0.762543 | ≈0 | 6372.048267 | — |
| 1 | 20 | 0° | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6377.503859 | — |
| 1 | 20 | 60° | 10 yr | Contact | — | — | Contact | n1 | 0.760780 | ≈0 | 6372.276790 | — |
| 1.2 | 15 | 0° | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6513.504239 | — |
| 1.2 | 15 | 60° | 10 yr | Contact | — | — | Contact | n1 | 0.763198 | ≈0 | 6373.542797 | — |
| 1.2 | 20 | 0° | 10 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6513.504237 | — |
| 1.2 | 20 | 60° | 10 yr | Contact | — | — | Contact | n1 | 0.762344 | ≈0 | 6373.640361 | — |
| 1 | 15 | 0° | 1 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6376.352496 | — |
| 1 | 15 | 60° | 1 yr | Contact | — | — | Contact | n1 | 0.762525 | ≈0 | 6369.959414 | — |
| 1 | 15 | 0° | Friction off | WD | — | — | WD | n1 | 0.001001 | 0.23082 | 6514.110202 | 2.22129 |
| 1 | 15 | 60° | Friction off | Contact | — | — | Contact | n1 | 0.762523 | ≈0 | 6372.005431 | — |
| 1 | 20 | 0° | Friction off | WD | — | — | WD | n1 | 0.001001 | 0.23070 | 6514.110202 | 2.22129 |
| 1 | 20 | 60° | Friction off | Contact | — | — | Contact | n1 | 0.760795 | ≈0 | 6372.461566 | — |
| 1 | 15 | 0° | 30 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6513.414788 | — |
| 1 | 15 | 60° | 30 yr | Contact | — | — | Contact | n1 | 0.762543 | ≈0 | 6372.040911 | — |
| 1 | 15 | 0° | 100 yr | Contact | — | — | Contact | n1 | 0.001000 | ≈0 | 6513.700823 | — |
| 1 | 15 | 60° | 100 yr | Contact | — | — | Contact | n1 | 0.762550 | ≈0 | 6372.020998 | — |
| 1 | 15 | 0° | Planet tides only | WD | WD | WD | WD | n2 | 0.001001 | 0.23082 | 6514.110202 | 2.22129 |
| 1 | 15 | 0° | 300 yr | WD | WD | WD | WD | n2 | 0.001000 | 0.17958 | 6514.110202 | 2.08798 |
| 1 | 15 | 60° | 300 yr | Contact | Contact | Contact | Contact | n2 | 0.762540 | ≈0 | 6372.009744 | — |
| 1 | 15 | 0° | 1,000 yr | WD | — | — | WD | n1 | 0.001000 | 0.22221 | 6514.110202 | 2.18696 |

</details>

### Data file

The CSV has 128 unique cases. It includes 208 numerical summary rows across the available resolutions and six additional outcome-only records from the refinement checks. The file includes sampled extrema, stop ages, contact-age differences and final orbital elements where the campaign reports them. The orbit histories and worker logs are stored separately on the run machine. Find it in **/data/planet_survival_results.csv**
