# Multi-Universal Gravity Simulation

A simulation for testing one idea: **universes are marbles resting on a shared sheet of space**, and the
gravity of neighbouring marbles, which lie beyond our horizon, reaches into ours.

The theory makes three claims:

1. **There is no dark matter.** The extra gravity we measure in galaxies comes from other universes.
2. **Universes collide.** Marbles drift together into dense wells, and those wells become new big bangs. The result looks like a neural network or a wavy sea of high and low points.
3. **Our edges are pulled unevenly**, because several neighbouring wells tug on different sides of our universe.

## Run it

Open `index.html` in a browser. It needs no build step. three.js loads from cdnjs.

- **Simulation 1: the sheet.** A flat sheet (switchable to a dented 3D view) with universes as marbles. Each marble's
  edge is its observable horizon. A *gravity reach* slider makes universes interact only when their edges are close.
  Mergers past a threshold trigger big bangs that throw matter back out as new universes. The teal ring shows the
  uneven (tidal) pull on our universe's edges.
- **Test 1: neighbours edge to edge.** A flat map of our universe (46.5 Gly radius) with neighbour universes that can
  touch our edge. Each universe's mass is spread over its whole disc. Click any galaxy to see its rotation curve
  compared with the flat measured curve. A gravity reach option adds a Yukawa fall-off.
- **Test 2: stacked sheets.** A hidden sheet a small distance away in an extra direction, holding a shadow galaxy under
  ours. "Find best fit" searches for the hidden mass and spread that best match the measured curve.

`analysis/tidal_check.py` (needs numpy) runs the original far-neighbour test in SI units and prints the numbers.

## What the tests show so far

| Claim | Result |
| --- | --- |
| Neighbours beside us explain dark matter | Fails, even edge to edge. The pull is 10⁵–10⁶× too weak. It averages to zero around a galaxy (Gauss's law). A reach limit turns the average outward. Gravity also travels at light speed, so matter we can't see can't pull on us either. |
| Stacked sheets explain dark matter | Can match the measured curve to within ~2 km/s, if the sheets are closer than ~10 kpc and the hidden galaxy is ~15–25× the visible mass and ~15× more spread out. That is a dark-matter halo on another sheet (brane-world "shadow matter"). |
| Collisions make big bangs | Plausible behaviour on the sheet. Real searches look for bubble-collision marks in the CMB (none confirmed). |
| Uneven pull on our edges | The testable part. It predicts direction-dependent expansion and aligned large-scale CMB anomalies. |
