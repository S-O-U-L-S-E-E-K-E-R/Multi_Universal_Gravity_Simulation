# Multi-Universal Gravity Simulation

A simulation for testing one idea: **universes are marbles resting on a shared sheet of space**, and the
gravity of neighbouring marbles, which lie beyond our horizon, reaches into ours.

The theory makes three claims:

1. **There is no dark matter.** The extra gravity we measure in galaxies comes from other universes.
2. **Universes collide.** Marbles drift together into dense wells, and those wells become new big bangs. The result looks like a neural network or a wavy sea of high and low points.
3. **Our edges are pulled unevenly**, because several neighbouring wells tug on different sides of our universe.

## Run it

Open `index.html` in a browser. It needs no build step. three.js loads from cdnjs.

- **Simulation 1: the sheet.** A 3D gravity sheet with marbles that merge. Mergers past a mass threshold trigger "big bangs".
  The teal ring around our universe shows the uneven (tidal) pull from the neighbours.
- **Test 1: the dark-matter claim.** A Milky Way–sized galaxy inside our universe, with 30 neighbour universes in random
  directions beyond our 46.5 Gly horizon. It compares the galaxy's rotation curve with the flat curve astronomers measure.

`analysis/tidal_check.py` (needs numpy) runs the same dark-matter test in SI units and prints the numbers.

## What the tests show so far

| Claim | Result |
| --- | --- |
| Other universes replace dark matter | Fails as stated. The pull is about 10⁵–10⁶× too weak. By Gauss's law it averages to zero around a galaxy. It is lopsided rather than centred on each galaxy. |
| Collisions make big bangs | Plausible behaviour on the sheet. Related real searches look for bubble-collision marks in the CMB (none confirmed). |
| Uneven pull on our edges | The testable part. It predicts direction-dependent expansion and aligned large-scale CMB anomalies. |

A version that could still explain dark matter: a neighbour universe that is close to ours in a *hidden* direction
(brane-world "shadow matter"). Its mass would then sit inside each galaxy's orbit.
