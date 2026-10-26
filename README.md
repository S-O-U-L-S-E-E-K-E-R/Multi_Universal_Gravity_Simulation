# Multi-Universal Gravity Simulation

A simulation for testing one idea: **universes are marbles resting on a shared sheet of space**, and the
gravity of neighbouring marbles, which lie beyond our horizon, reaches into ours.

The theory makes three claims:

1. **There is no dark matter.** The extra gravity we measure in galaxies comes from other universes.
2. **Universes collide.** Marbles drift together into dense wells, and those wells become new big bangs. The result looks like a neural network or a wavy sea of high and low points.
3. **Our edges are pulled unevenly**, because several neighbouring wells tug on different sides of our universe.

## Where the theory stands

| Claim | Status |
| --- | --- |
| Another universe's matter, mixed into our space and interacting only by gravity, is dark matter | Works: passes rotation curves (model and 149 real galaxies), clusters and the Bullet Cluster |
| Our universe falling into a bigger well as a whole | Allowed at any strength, but undetectable from inside |
| Neighbours' wells (from outside) as dark matter | Ruled out: averages to zero around a galaxy and is far too weak |
| The stretch from a bigger well as the accelerating expansion | Ruled out by supernovae and by the CMB (needed stretch is 10⁴–10⁵× the allowed one) |
| Big bangs at extreme density | Plausible (bounce / black-hole cosmology), untested |
| Uneven edges, region-dependent cores | Hints only, not significant |

**In what world it works:** our universe ignited inside a well of another universe's matter, which still fills our
space as dark matter (about 5× ours, gravity-only, non-cooling); dark energy is the stretching of the sheet that
carries all universes, not a slope; neighbours are far or light (one as heavy as ours more than ~19 of our radii
away); and big bangs need extreme density to bounce instead of forming black holes.

## Run it

Open `index.html` in a browser. It needs no build step. three.js loads from cdnjs.

- **Simulation 1: the sheet.** Each universe is a collection of freely moving weights (dense centre, thin edge) that
  starts with its own big bang, so its edge grows and its centre thins over time. Every weight pulls on every other
  (particle-mesh gravity), so neighbours can dent and stretch a universe's shape. Big bangs are set off by density
  alone: any patch denser than the chosen threshold (measured against a newborn universe's core) ignites. A universe
  whose own core falls back in and re-ignites keeps its identity. Total mass is conserved.
- **Well profile.** A side view of one universe's well as it expands (fixed mass, centrally concentrated). The centre
  gets shallower and the peak pull weakens fast, while the pull at a distant point barely changes until the
  universe's matter spreads out to reach it.
- **Test 1: neighbours edge to edge.** A flat map of our universe (46.5 Gly radius) with neighbour universes that can
  touch our edge. Each universe's mass is spread over its whole disc. Click any galaxy to see its rotation curve
  compared with the flat measured curve. A gravity reach option adds a Yukawa fall-off.
- **Test 1b: overlapping neighbour.** An expanded neighbour whose thin matter fills our galaxy. Its density is set
  relative to our universe's average (about 4 × 10⁻²⁸ kg/m³).
- **Simulation 2: galaxies forming with a neighbour's matter.** Our matter and a neighbour universe's matter share
  expanding space, each starting with its own independent random lumps, coupled only by gravity. It measures whether the sheets line up, whether
  hidden halos form around our galaxies, and how flat the rotation curves come out. It includes a control run with no
  gravity across the gap. `twosheet.js` logic is inlined into `index.html`.
- **Test 2: a neighbour's matter inside our galaxy.** A cloud of another universe's matter gathered around our galaxy,
  interacting only through gravity, optionally kept out of the very centre (a core). "Find best fit" searches for the hidden mass and spread that best match the measured curve.

`analysis/tidal_check.py` (needs numpy) runs the original far-neighbour test in SI units and prints the numbers.
`analysis/anomaly_axes.py` (needs numpy) tests whether published cosmic anomaly directions share one axis.
`analysis/dark_matter_tests.py` (needs numpy) is the test suite: seven versions of "other universes do dark matter's
job", each fitted on galaxy rotation and then scored, without refitting, against the speed–mass rule, every-galaxy
coverage, galaxy-cluster mass and the Bullet Cluster.

`analysis/sparc_tests.py` and `analysis/sparc_per_galaxy_gap.py` test the versions against 149 real galaxies from
the SPARC catalogue (download links in the script headers).

## Real-galaxy result (SPARC)

With one set of settings for all galaxies, neighbour matter in our wells cuts the rotation-curve miss from 0.57 dex (visible
matter only) to 0.18 dex, and keeping it out of galaxy centres (a ~3 kpc core) fits clearly better than none. Fitted galaxy by
galaxy, though, no single core size fits most of them (48% vs the ~68% a universal gap needs): core sizes grow weakly
with galaxy size. Star-explosion gas loss added to the two-sheet simulation did not produce the hidden-matter
trend this version needs. Core sizes that vary by place predict that galaxies near each other share core sizes. Tested with
`analysis/wavy_gap_test.py` (positions in `data/sparc_coords.json`, from the CDS Sesame resolver): nearby galaxies'
core sizes are ~25% more similar than random pairs (p = 0.003), but after controlling for distance from us the same
effect is only p = 0.11. Suggestive, not established.

## Accelerating expansion: a tide from a bigger well?

`analysis/supernova_tide_test.py` fits 1,309 Pantheon+ supernovae. A tide (stretching along one axis, squeezing at
right angles, no dark energy on average) does better than no acceleration at all but is decisively worse than dark
energy that is the same in every direction (BIC +68.5 with 30% matter; +44 even with no matter). Adding a tide on
top of dark energy is not supported either (BIC +18).

`analysis/supernova_motion_rotation_test.py` adds a lopsided (near-side stronger) pull, a rotating universe, and our
own motion left free (raw Sun-centred redshifts). Dark energy stays best in both frames (lopsided well alone +29/+30,
rotation alone +18/+22, rotation + lopsided well +24 BIC). The supernovae recover our own motion independently
(250–350 km/s, 19–59° from the CMB direction). CMB limits on uneven expansion and rotation (Saadeh et al. 2016) are
far stricter still; falling into a bigger well as a whole remains undetectable from inside.

## Simulation 3: our universe on a slope

The page's third simulation puts our universe on the slope of one or two neighbours' wells (mass, distance, angle).
From outside the whole universe slides toward the well; from inside only the uneven stretch shows. Turned around, the
CMB caps how heavy a neighbour can be at each distance (`analysis/neighbour_limits.py`): a neighbour as heavy as our
universe must be more than ~19 of our radii away, and at every distance the mass needed to replace dark energy is
~10⁴–10⁵ times what the CMB allows.

## Test suite result

Only a neighbour universe's matter that has expanded into our space and gathered in our galaxies' wells passes all five tests, and only if small galaxies
have also lost part of their ordinary matter (a process the simulations here don't model). Gravity leaking into the
sheet below a fixed acceleration (MOND-like) matches galaxies but fails clusters and the Bullet Cluster. Neighbours
beside us, sitting in a bigger well, and overlapping neighbours fail every test.

## What the tests show so far

| Claim | Result |
| --- | --- |
| Neighbours beside us explain dark matter | Fails, even edge to edge. The pull is 10⁵–10⁶× too weak. It averages to zero around a galaxy (Gauss's law). A reach limit turns the average outward. (Light-speed limits on gravity don't rule it out if the neighbours' wells predate our universe.) |
| An expanded neighbour overlapping us explains it | Fails. It needs to be ~2 × 10⁵ times denser than our universe, and even then the curve rises instead of staying flat (best fit off by ~42 km/s). |
| A neighbour's matter in our space explains dark matter | Can match the measured curve to within ~2 km/s if the cloud of neighbour matter around a galaxy is ~15–25× its visible mass and ~15× more spread out: a dark-matter halo made of another universe's matter. |
| Neighbour matter gathers in our wells on its own | Yes, in Simulation 2: the two kinds of matter line up (~0.93–0.95 vs ~0.02 in the control), neighbour clouds are 5–10× wider than galaxies, and curves come out flatter (0.60–0.72 vs ~0.47). Without the neighbour's matter, no galaxies form at all in that time. |
| Extreme density makes big bangs | Plausible behaviour on the sheet. Real searches look for bubble-collision marks in the CMB (none confirmed). |
| Uneven pull on our edges | Testable, but not supported yet. Six published anomaly axes are no more aligned than random (p ≈ 0.5). The bulk flow, Hubble-rate anisotropy and dark flow sit within ~20° of (l, b) ≈ (289°, 3°), but that grouping is consistent with chance and lies near the Great Attractor. |
