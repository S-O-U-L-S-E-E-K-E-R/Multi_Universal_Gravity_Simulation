"""Can gravity from other universes reproduce what dark matter does?

Each version of the idea gets its best-fit global settings from galaxy rotation (the radial acceleration relation),
then has to pass the other tests with those same settings, without refitting.

Targets (all measured):
  * Radial acceleration relation (RAR): g_obs = g_bar / (1 - exp(-sqrt(g_bar / 1.2e-10))), observed scatter ~0.13 dex
    (McGaugh, Lelli & Schombert 2016, 153 SPARC galaxies).
  * Baryonic Tully-Fisher relation (BTFR): M_b proportional to v_flat^slope with slope ~3.85 (Lelli et al. 2016).
  * Galaxy clusters: total mass ~6x the visible (baryonic) mass inside ~1 Mpc (cosmic baryon fraction ~0.16).
  * Bullet Cluster: lensing mass peaks sit with the galaxies, ~0.7 Mpc from the colliding gas (Clowe et al. 2006).

Usage: python3 analysis/dark_matter_tests.py [--json out.json]
"""
import json
import sys

import numpy as np

G = 6.674e-11
MSUN = 1.989e30
KPC = 3.086e19
MPC = 1e3 * KPC
G_DAGGER = 1.2e-10  # RAR acceleration scale, m/s^2
RAR_SCATTER = 0.13  # observed scatter, dex
rng = np.random.default_rng(2026)


# ---------------------------------------------------------------- galaxies
def make_galaxies(n=60):
    """Model disc galaxies spanning dwarfs to giants. Size grows with mass (roughly as observed)."""
    mb = 10 ** rng.uniform(8.0, 11.5, n) * MSUN
    rd = 3.0 * KPC * (mb / (5e10 * MSUN)) ** 0.33 * 10 ** rng.normal(0, 0.1, n)
    normals = rng.normal(size=(n, 3))
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return [dict(mb=mb[i], rd=rd[i], n=normals[i]) for i in range(n)]


def radii(g):
    return np.linspace(0.5, 8, 24) * g["rd"]


def g_bar(g, r):
    """Pull from visible matter (exponential disc, enclosed-mass approximation)."""
    x = r / g["rd"]
    m = g["mb"] * (1 - (1 + x) * np.exp(-x))
    return G * m / r**2


def g_target(gb):
    return gb / (1 - np.exp(-np.sqrt(gb / G_DAGGER)))


# ---------------------------------------------------------------- the versions of the idea
# Each takes (galaxy, radii, params) and returns the total inward pull. Params are global: one universe, one set.
def tidal_tensor(amplitude):
    """A neighbour universe's tide: traceless, stretching along x. Every galaxy in our universe feels the same one."""
    return amplitude * np.diag([2.0, -1.0, -1.0])


def mech_A(g, r, p):  # neighbours beside us: only their tide reaches a galaxy's orbits
    T = tidal_tensor(p["tide"])
    n = g["n"]
    # average inward pull around a ring in the galaxy's disc: +r * (n.T.n) / 2 (the trace is zero)
    extra = r * (n @ T @ n) / 2
    return g_bar(g, r) + extra


def mech_B(g, r, p):  # sitting inside a big outside well: a uniform pull moves the whole galaxy together
    return g_bar(g, r) + 0.0


def mech_C(g, r, p):  # an expanded neighbour overlapping us at even density rho
    return g_bar(g, r) + 4 / 3 * np.pi * G * p["rho"] * r


SIM2_RATIO = 5.0  # hidden-to-visible mass ratio measured in Simulation 2 (71 galaxies, no trend with mass)


def mech_D(g, r, p):  # stacked hidden sheet: a hidden galaxy (mass f * visible, spread k * disc size) a gap h away
    f = p.get("f", SIM2_RATIO)
    if "alpha" in p:  # small galaxies have lost part of their visible matter, so they hold relatively more hidden matter
        f = p["f0"] * (g["mb"] / (1e10 * MSUN)) ** (-p["alpha"])
    mh, s, h = f * g["mb"], p["k"] * g["rd"], p["h"] * KPC
    rho_ = np.sqrt(r**2 + h**2)
    return g_bar(g, r) + G * mh / (rho_ + s) ** 2 * (r / rho_)


def mech_E1(g, r, p):  # gravity leaks into the sheet beyond a fixed distance: falls as 1/r instead of 1/r^2 out there
    gb = g_bar(g, r)
    rc = p["rc"] * KPC
    return np.where(r < rc, gb, gb * r / rc)


def mech_E2(g, r, p):  # gravity leaks into the sheet below a fixed acceleration a0 (MOND-like)
    gb = g_bar(g, r)
    return gb * (0.5 + np.sqrt(0.25 + p["a0"] / gb))


MECHS = {
    "A": ("Neighbour universes beside us (their tide)", mech_A, {"tide": np.concatenate([-np.logspace(-40, -28, 60), np.logspace(-40, -28, 60)])}),
    "B": ("We sit inside a big outside well", mech_B, {}),
    "C": ("An expanded neighbour overlapping us", mech_C, {"rho": np.logspace(-30, -18, 121)}),
    "Da": ("A stacked hidden sheet, as Simulation 2 forms it", mech_D, {"k": np.linspace(1, 40, 79), "h": [0, 2, 5, 10, 30]}),
    "Db": ("A stacked hidden sheet + small galaxies lose visible matter", mech_D, {"f0": np.linspace(2, 30, 29), "alpha": np.linspace(0, 0.8, 17), "k": np.linspace(2, 30, 15), "h": [0.0]}),
    "E1": ("Gravity leaks into the sheet beyond a fixed distance", mech_E1, {"rc": np.logspace(-0.5, 2.5, 61)}),
    "E2": ("Gravity leaks into the sheet below a fixed acceleration", mech_E2, {"a0": np.logspace(-12, -8, 81)}),
}


def rar_rms(fn, gals, p):
    res = []
    for g in gals:
        r = radii(g)
        pred = fn(g, r, p)
        tgt = g_target(g_bar(g, r))
        bad = pred <= 0  # pushed outward: no orbit possible; count as a large miss
        lr = np.log10(np.where(bad, 1e-3 * tgt, pred) / tgt)
        res.append(lr)
    res = np.concatenate(res)
    return float(np.sqrt(np.mean(res**2)))


def fit(fn, gals, grid):
    keys = list(grid)
    if not keys:
        return {}, rar_rms(fn, gals, {})
    best = (np.inf, None)
    for combo in np.array(np.meshgrid(*[grid[k] for k in keys], indexing="ij")).reshape(len(keys), -1).T:
        p = dict(zip(keys, combo))
        e = rar_rms(fn, gals, p)
        if e < best[0]:
            best = (e, p)
    return best[1], best[0]


def btfr_slope(fn, gals, p):
    """Slope of log(visible mass) against log(outer orbital speed)."""
    lm, lv = [], []
    for g in gals:
        r = np.array([6 * g["rd"]])
        a = fn(g, r, p)[0]
        if a <= 0:
            continue
        lm.append(np.log10(g["mb"] / MSUN))
        lv.append(np.log10(np.sqrt(a * r[0]) / 1e3))
    if len(lm) < 5:
        return float("nan"), 0
    return float(np.polyfit(lv, lm, 1)[0]), len(lm)


def per_galaxy_pass(fn, gals, p):
    ok = 0
    for g in gals:
        r = radii(g)
        pred = fn(g, r, p)
        tgt = g_target(g_bar(g, r))
        if np.all(pred > 0) and np.sqrt(np.mean(np.log10(pred / tgt) ** 2)) < RAR_SCATTER:
            ok += 1
    return ok / len(gals)


def flatness(fn, gals, p):
    vals = []
    for g in gals:
        r = np.array([3 * g["rd"], 8 * g["rd"]])
        a = fn(g, r, p)
        if np.all(a > 0):
            v = np.sqrt(a * r)
            vals.append(v[1] / v[0])
    return float(np.median(vals)) if vals else float("nan")


# ---------------------------------------------------------------- galaxy clusters
def cluster_ratio(key, p):
    """Mass discrepancy (dynamical / visible) inside 1 Mpc for a 1e14 Msun gas cluster, using galaxy-fit settings."""
    r = 1.0 * MPC
    rc, mgas, mstar = 0.15 * MPC, 1.2e14 * MSUN, 1.5e13 * MSUN
    x = np.linspace(1e-3, 1, 4000) * r  # beta = 2/3 gas profile, normalised to mgas inside 1 Mpc
    dens = 1 / (1 + (x / rc) ** 2)
    mb = mgas + mstar
    gb = G * mb / r**2
    if key == "A":
        return 1.0  # a tide averages to zero over the cluster: no net extra mass
    if key == "B":
        return 1.0
    if key == "C":
        return 1 + (4 / 3 * np.pi * G * p["rho"] * r) / gb
    if key in ("Da", "Db"):
        # clusters are too massive to lose their gas, so they keep the ratio the two sheets formed with (Simulation 2)
        return 1 + SIM2_RATIO
    if key == "E1":
        rcut = p["rc"] * KPC
        return r / rcut if r > rcut else 1.0
    if key == "E2":
        return 0.5 + np.sqrt(0.25 + p["a0"] / gb)
    return float("nan")


# ---------------------------------------------------------------- Bullet Cluster
def bullet_offset(key):
    """Where does the lensing mass of one sub-cluster peak, measured from its gas (Mpc)?
    After the collision the gas (about 85% of the visible mass) is left behind near the middle and the galaxies
    (15%) have moved on, 0.7 Mpc ahead. Observed: the lensing peak sits with the galaxies."""
    gas, stars, sep = 0.85, 0.15, 0.7
    if key in ("Da", "Db"):
        hidden = SIM2_RATIO  # collisionless hidden matter passes through with the galaxies
        return (stars + hidden) * sep / (gas + stars + hidden)
    # every other version only adds pull in proportion to (or independent of) visible matter, so the
    # lensing peak follows the visible mass, which is dominated by the gas
    return stars * sep / (gas + stars)


def main():
    gals = make_galaxies()
    results = {}
    for key, (label, fn, grid) in MECHS.items():
        p, rms = fit(fn, gals, grid)
        slope, n_ok = btfr_slope(fn, gals, p)
        res = dict(
            label=label,
            params={k: float(v) for k, v in p.items()},
            rar_rms_dex=rms,
            btfr_slope=slope,
            galaxies_fit_fraction=per_galaxy_pass(fn, gals, p),
            flatness=flatness(fn, gals, p),
            cluster_ratio=float(cluster_ratio(key, p)),
            bullet_offset_mpc=bullet_offset(key),
        )
        res["passes"] = dict(
            rotation=rms < RAR_SCATTER,
            btfr=bool(np.isfinite(slope) and 3.5 <= slope <= 4.3),
            every_galaxy=res["galaxies_fit_fraction"] >= 0.9,
            clusters=4.5 <= res["cluster_ratio"] <= 8,
            bullet=res["bullet_offset_mpc"] >= 0.4,
        )
        # RAR points for plotting: one curve per galaxy, coarse
        pts = []
        for g in gals[::3]:
            r = radii(g)[::2]
            gb = g_bar(g, r)
            pr = fn(g, r, p)
            pts.append([[float(np.log10(a)), float(np.log10(b)) if b > 0 else None] for a, b in zip(gb, pr)])
        res["rar_points"] = pts
        results[key] = res

    ref = []
    for lg in np.linspace(-12.5, -8.5, 41):
        ref.append([float(lg), float(np.log10(g_target(10**lg)))])

    print("Can other-universe gravity reproduce dark matter?  (galaxy-fit settings reused for every test)\n")
    hdr = f"{'':4s}{'version':55s}{'RAR miss':>9s}{'BTFR':>7s}{'gal ok':>8s}{'flat':>6s}{'cluster':>9s}{'bullet':>8s}"
    print(hdr)
    print(f"{'':4s}{'target':55s}{'<0.13':>9s}{'~3.85':>7s}{'>=90%':>8s}{'~1':>6s}{'~6x':>9s}{'~0.7':>8s}")
    for key, r in results.items():
        ps = r["passes"]
        mark = lambda b: "✓" if b else "✗"
        print(f"{key:4s}{r['label']:55s}{r['rar_rms_dex']:7.2f}{mark(ps['rotation']):>2s}"
              f"{r['btfr_slope']:6.2f}{mark(ps['btfr'])}{100*r['galaxies_fit_fraction']:6.0f}%{mark(ps['every_galaxy'])}"
              f"{r['flatness']:6.2f}{r['cluster_ratio']:8.1f}{mark(ps['clusters'])}{r['bullet_offset_mpc']:7.2f}{mark(ps['bullet'])}")
        print(f"{'':4s}  best settings: {r['params']}")
    print(f"\nTargets: RAR miss = rms gap from the measured relation (dex); BTFR = slope of log M vs log v;")
    print("gal ok = share of galaxies matched within 0.13 dex; flat = v(8 Rd) / v(3 Rd); cluster = mass discrepancy at 1 Mpc;")
    print("bullet = lensing peak offset from the gas (Mpc).")

    if "--json" in sys.argv:
        out = sys.argv[sys.argv.index("--json") + 1]
        with open(out, "w") as f:
            json.dump({"results": results, "rar_reference": ref}, f)
        print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
