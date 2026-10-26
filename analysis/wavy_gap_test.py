"""Wavy-gap test: if the hidden sheet ripples, the gap changes from place to place, so galaxies that are close
together in space should have similar gaps (core sizes), whatever their own size.

Inputs: per-galaxy gap fits (from sparc_per_galaxy_gap.py), sky coordinates (from the CDS Sesame resolver) and the
SPARC distances. The known trend (bigger galaxies prefer bigger gaps) is removed first, so it can't fake a signal.
Significance comes from shuffling the gaps among the galaxies (positions fixed), and the "best of several distance
scales" is corrected by taking the same best-of-scales in every shuffle.

Usage: python3 analysis/wavy_gap_test.py <pergal.json> <coords.json> <SPARC_Lelli2016c.mrt>
"""
import json
import sys

import numpy as np

rng = np.random.default_rng(12)
pergal, coords, mrt = sys.argv[1:4]
dist = {}
for line in open(mrt):
    t = line.split()
    if len(t) >= 18 and t[1].isdigit():
        dist[t[0]] = (float(t[2]), float(t[3]))  # distance and its error, Mpc
rows = [x for x in json.load(open(pergal)) if x["hi"] - x["lo"] <= 3 and x["hi"] < 13]  # gap pinned down
cd = json.load(open(coords))
rows = [x for x in rows if cd.get(x["name"]) and x["name"] in dist]

ra = np.radians([cd[x["name"]][0] for x in rows]); dec = np.radians([cd[x["name"]][1] for x in rows])
d = np.array([dist[x["name"]][0] for x in rows])
xyz = np.stack([d * np.cos(dec) * np.cos(ra), d * np.cos(dec) * np.sin(ra), d * np.sin(dec)], 1)
h = np.array([x["h"] for x in rows]); lrd = np.log10([x["rd"] for x in rows])
ld = np.log10(d)
print(f"gap vs distance from us: correlation {np.corrcoef(ld, h)[0,1]:.2f}, slope {np.polyfit(ld, h, 1)[0]:.2f} kpc per factor 10 in distance")
CONTROL_DISTANCE = "--control-distance" in sys.argv
if CONTROL_DISTANCE:  # remove both trends: galaxy size and distance from us (a measurement effect)
    A = np.stack([np.ones(len(h)), lrd, ld], 1)
    resid = h - A @ np.linalg.lstsq(A, h, rcond=None)[0]
else:
    resid = h - np.polyval(np.polyfit(lrd, h, 1), lrd)  # gap left over after the size trend
DIST_BINS = np.digitize(ld, np.quantile(ld, [0.2, 0.4, 0.6, 0.8]))  # five equal-count distance groups
n = len(rows)
sep3 = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
cosang = np.clip(np.sin(dec)[:, None] * np.sin(dec)[None] + np.cos(dec)[:, None] * np.cos(dec)[None] * np.cos(ra[:, None] - ra[None]), -1, 1)
sepang = np.degrees(np.arccos(cosang))
iu = np.triu_indices(n, 1)


def run(sep, scales, label, unit):
    s = sep[iu]
    masks = [s < sc for sc in scales]
    def stats(r):
        dr = np.abs(r[:, None] - r[None])[iu]
        return np.array([dr[m].mean() / dr.mean() if m.sum() >= 5 else np.nan for m in masks])
    obs = stats(resid)
    def shuffle(r):
        if not CONTROL_DISTANCE:
            return rng.permutation(r)
        out = r.copy()  # shuffle only within groups of similar distance, so any distance bias is kept in every shuffle
        for b in np.unique(DIST_BINS):
            i = np.where(DIST_BINS == b)[0]; out[i] = r[rng.permutation(i)]
        return out
    sims = np.array([stats(shuffle(resid)) for _ in range(20000)])
    print(f"\n{label}: {n} galaxies with a pinned-down gap")
    pvals = []
    for i, sc in enumerate(scales):
        p = np.mean(sims[:, i] <= obs[i]); pvals.append(p)
        print(f"  pairs closer than {sc:5g} {unit}: {masks[i].sum():4d} pairs, gap difference {obs[i]:.2f}x the average  (chance p = {p:.3f})")
    # best-of-scales, corrected: how often does a shuffle reach a p this small at ANY scale?
    sim_p = np.array([[np.mean(sims[:, i] <= sims[j, i]) for i in range(len(scales))] for j in range(0, 20000, 10)])
    pmin = np.nanmin(pvals); corrected = np.mean(np.nanmin(sim_p, 1) <= pmin)
    print(f"  smallest p {pmin:.3f}; corrected for trying {len(scales)} scales: p = {corrected:.3f}")
    return dict(scales=list(scales), ratio=obs.tolist(), p=pvals, p_corrected=float(corrected),
                pairs=[int(m.sum()) for m in masks])


res = dict(n=n,
           space=run(sep3, [1, 2, 3, 5, 8, 12, 20], "Distance in space (Mpc)", "Mpc"),
           sky=run(sepang, [2, 5, 10, 20, 40], "Angle on the sky (deg)", "deg"))
print("\n(ratio < 1 means nearby galaxies have more similar gaps than random pairs; a wavy gap predicts < 1)")
outs = [a for a in sys.argv[4:] if not a.startswith("--")]
if outs:
    json.dump(res, open(outs[0], "w"))
