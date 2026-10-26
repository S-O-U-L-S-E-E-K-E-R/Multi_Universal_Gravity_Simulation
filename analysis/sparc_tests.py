"""Test versions of the theory against real galaxies: the SPARC catalogue (Lelli, McGaugh & Schombert 2016).

Every version gets ONE set of global settings for all galaxies, so it has to predict each galaxy's measured rotation
curve from that galaxy's visible matter alone.

Data (not included here; public):
  https://astroweb.case.edu/SPARC/Rotmod_LTG.zip        -> unzip into <dir>/rotmod/
  https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt  -> <dir>/SPARC_Lelli2016c.mrt

Usage: python3 analysis/sparc_tests.py <dir> [--json out.json]
"""
import glob
import json
import os
import sys
from itertools import product

import numpy as np

G = 4.30091e-6  # kpc (km/s)^2 / Msun
UPS_DISK, UPS_BUL = 0.5, 0.7  # standard 3.6-micron mass-to-light ratios


def load(d):
    meta = {}
    for line in open(os.path.join(d, "SPARC_Lelli2016c.mrt")):
        t = line.split()  # name T D e_D f_D Inc e_Inc L e_L Reff SBeff Rdisk SBdisk MHI RHI Vflat e_Vflat Q Ref
        if len(t) < 18 or not t[1].isdigit():
            continue
        try:
            meta[t[0]] = dict(inc=float(t[5]), L=float(t[7]) * 1e9, rd=float(t[11]), mhi=float(t[13]) * 1e9, q=int(t[17]))
        except ValueError:
            continue
    gals = []
    for f in sorted(glob.glob(os.path.join(d, "rotmod", "*_rotmod.dat"))):
        name = os.path.basename(f).replace("_rotmod.dat", "")
        m = meta.get(name)
        if not m or m["q"] >= 3 or m["inc"] < 30:
            continue
        a = np.loadtxt(f, comments="#")
        r, vobs, ev, vg, vd, vb = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
        vbar2 = vg * np.abs(vg) + UPS_DISK * vd * np.abs(vd) + UPS_BUL * vb * np.abs(vb)
        ok = (r > 0) & (vobs > 0) & (vbar2 > 0)
        if ok.sum() < 5:
            continue
        mb = UPS_DISK * m["L"] + 1.33 * m["mhi"]
        gals.append(dict(name=name, r=r[ok], gobs=vobs[ok] ** 2 / r[ok], gbar=vbar2[ok] / r[ok],
                         sig=np.maximum(2 * ev[ok] / (vobs[ok] * np.log(10)), 0.02), mb=mb, rd=max(m["rd"], 0.1)))
    return gals


# ---- versions of the theory (all accelerations in (km/s)^2/kpc) ----
A0_UNIT = 3.0857e19 / 1e6  # converts m/s^2 to (km/s)^2/kpc


def newton(g, p):
    return g["gbar"]


def leak_accel(g, p):  # E2: gravity leaks into the sheet below a fixed acceleration
    a0 = p["a0"] * A0_UNIT
    return g["gbar"] * (0.5 + np.sqrt(0.25 + a0 / g["gbar"]))


def stacked(g, p):  # D-b: hidden galaxy on a stacked sheet, mass f(M) * visible, spread k * disc size, gap h
    f = p["f0"] * (g["mb"] / 1e10) ** (-p["alpha"])
    mh, s, h, r = f * g["mb"], p["k"] * g["rd"], p["h"], g["r"]
    rho = np.sqrt(r**2 + h**2)
    return g["gbar"] + G * mh / (rho + s) ** 2 * (r / rho)


def score(fn, gals, p):
    """Each galaxy counts equally: mean over galaxies of the mean squared log miss, weighted by measurement error."""
    tot, per = 0.0, []
    for g in gals:
        pred = fn(g, p)
        d = np.log10(pred / g["gobs"])
        w = 1 / (g["sig"] ** 2 + 0.05**2)  # 0.05 dex floor for distance/inclination systematics
        per.append(float(np.sqrt(np.sum(w * d**2) / np.sum(w))))
        tot += np.sum(w * d**2) / np.sum(w)
    return tot / len(gals), per


def fit(fn, gals, grid):
    keys = list(grid)
    best = (np.inf, None)
    for combo in product(*[grid[k] for k in keys]):
        p = dict(zip(keys, combo))
        s, _ = score(fn, gals, p)
        if s < best[0]:
            best = (s, p)
    return best


def summarise(label, fn, gals, p):
    s, per = score(fn, gals, p)
    per = np.array(per)
    return dict(label=label, params=p, rms_dex=float(np.sqrt(s)), median_galaxy_dex=float(np.median(per)),
                good_fraction=float(np.mean(per < 0.1)), per_galaxy=per.tolist())


def main():
    d = sys.argv[1]
    gals = load(d)
    print(f"SPARC galaxies used: {len(gals)}  (quality 1-2, inclination >= 30 deg), {sum(len(g['r']) for g in gals)} points\n")
    out = {}
    out["N"] = summarise("Visible matter only", newton, gals, {})
    s, p = fit(leak_accel, gals, {"a0": np.linspace(0.6e-10, 2.0e-10, 57)})
    out["E2"] = summarise("Gravity leaks into the sheet below a fixed acceleration", leak_accel, gals, p)
    coarse = {"f0": [4, 6, 8, 10, 12, 15, 18, 22, 27, 33], "alpha": [0, 0.1, 0.2, 0.3, 0.4, 0.5],
              "k": [2, 3, 4, 6, 8, 11, 15, 20], "h": [0.0]}
    s, p = fit(stacked, gals, coarse)
    fine = {"f0": np.linspace(p["f0"] * 0.7, p["f0"] * 1.3, 9), "alpha": np.linspace(max(0, p["alpha"] - 0.1), p["alpha"] + 0.1, 9),
            "k": np.linspace(p["k"] * 0.7, p["k"] * 1.3, 9), "h": [0.0]}
    s, p0 = fit(stacked, gals, fine)
    out["Db"] = summarise("Stacked hidden sheet, no gap (hidden pull peaks at each centre)", stacked, gals, p0)
    # the gap: a stacked sheet a distance h away softens the hidden pull inside ~h, giving every galaxy a core
    scan = []
    for h in [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8]:
        g2 = {"f0": np.linspace(p0["f0"] * 0.6, p0["f0"] * 1.8, 13), "alpha": np.linspace(max(0, p0["alpha"] - 0.15), p0["alpha"] + 0.15, 7),
              "k": np.linspace(p0["k"] * 0.5, p0["k"] * 1.5, 11), "h": [h]}
        sh, ph = fit(stacked, gals, g2)
        scan.append(dict(h=h, rms_dex=float(np.sqrt(sh)), params=ph))
        print(f"  gap {h:5.2f} kpc: miss {np.sqrt(sh):.4f} dex  (f0={ph['f0']:.1f}, alpha={ph['alpha']:.2f}, k={ph['k']:.1f})")
    best = min(scan, key=lambda x: x["rms_dex"])
    out["Dgap"] = summarise(f"Stacked hidden sheet, best gap {best['h']} kpc", stacked, gals, best["params"])
    print()
    for key, r in out.items():
        print(f"{key:5s} {r['label']:66s} miss {r['rms_dex']:.3f} dex | median galaxy {r['median_galaxy_dex']:.3f} | "
              f"galaxies within 0.1 dex: {100*r['good_fraction']:.0f}%  {({k: round(float(v), 3) for k, v in r['params'].items()})}")
    if "--json" in sys.argv:
        with open(sys.argv[sys.argv.index("--json") + 1], "w") as f:
            json.dump(dict(results=out, gap_scan=scan, n_gal=len(gals),
                           names=[g["name"] for g in gals], mb=[g["mb"] for g in gals]), f, default=float)


if __name__ == "__main__":
    main()
