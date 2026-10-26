# Multi-Universal Gravity Simulation

A simulation for testing one idea: **universes are marbles resting on a shared sheet of space**, and the
gravity of neighbouring marbles, which lie beyond our horizon, reaches into ours.

The theory makes three claims:

1. **There is no dark matter.** The extra gravity we measure in galaxies comes from other universes.
2. **Universes collide.** Marbles drift together into dense wells, and those wells become new big bangs. The result looks like a neural network or a wavy sea of high and low points.
3. **Our edges are pulled unevenly**, because several neighbouring wells tug on different sides of our universe.

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
- **Simulation 2: two sheets forming galaxies.** Two stacked sheets in expanding space, each starting with its own
  independent random lumps, coupled only by gravity across the gap. It measures whether the sheets line up, whether
  hidden halos form around our galaxies, and how flat the rotation curves come out. It includes a control run with no
  gravity across the gap. `twosheet.js` logic is inlined into `index.html`.
- **Test 2: stacked sheets.** A hidden sheet a small distance away in an extra direction, holding a shadow galaxy under
  ours. "Find best fit" searches for the hidden mass and spread that best match the measured curve.

`analysis/tidal_check.py` (needs numpy) runs the original far-neighbour test in SI units and prints the numbers.

## What the tests show so far

| Claim | Result |
| --- | --- |
| Neighbours beside us explain dark matter | Fails, even edge to edge. The pull is 10⁵–10⁶× too weak. It averages to zero around a galaxy (Gauss's law). A reach limit turns the average outward. (Light-speed limits on gravity don't rule it out if the neighbours' wells predate our universe.) |
| An expanded neighbour overlapping us explains it | Fails. It needs to be ~2 × 10⁵ times denser than our universe, and even then the curve rises instead of staying flat (best fit off by ~42 km/s). |
| Stacked sheets explain dark matter | Can match the measured curve to within ~2 km/s, if the sheets are closer than ~10 kpc and the hidden galaxy is ~15–25× the visible mass and ~15× more spread out. That is a dark-matter halo on another sheet (brane-world "shadow matter"). |
| Stacked sheets form on their own | Yes, in Simulation 2: the sheets line up (~0.9 vs ~0.02 in the control), hidden halos are several times wider than galaxies, and curves come out flatter (~0.8 vs ~0.47). Without the hidden sheet, no galaxies form at all in that time. |
| Extreme density makes big bangs | Plausible behaviour on the sheet. Real searches look for bubble-collision marks in the CMB (none confirmed). |
| Uneven pull on our edges | The testable part. It predicts direction-dependent expansion and aligned large-scale CMB anomalies. |
