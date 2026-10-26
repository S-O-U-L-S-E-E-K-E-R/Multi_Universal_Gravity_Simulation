"""Is the accelerating expansion really a tide from a bigger well (our universe stretched like spaghetti)?

A tide stretches along one axis and squeezes at right angles. Averaged over the whole sky it adds nothing (a tide's
trace is zero). So "tide instead of dark energy" predicts acceleration along one axis and extra slowing at right
angles, while dark energy predicts the same acceleration everywhere. Type Ia supernovae across the sky decide.

Each direction n gets its own expansion history with an effective dark-energy term
    Omega_L(n) = Omega_L0 + eps * P2(n . axis),   P2(x) = (3x^2 - 1) / 2
(eps > 0: extra acceleration along the axis, slowing at right angles), and curvature 1 - Omega_m - Omega_L(n).

Models (absolute magnitude / H0 offset fitted analytically in all):
  matter only    : Omega_L0 = 0, eps = 0, Omega_m free
  dark energy    : flat, Omega_L0 = 1 - Omega_m, eps = 0
  tide only      : Omega_L0 = 0, eps and axis free  (the theory: no dark energy, only stretching)
  dark energy + tide

Data: Pantheon+ (Scolnic et al. 2022; Brout et al. 2022), diagonal errors only.
  https://github.com/PantheonPlusSH0ES/DataRelease  ->  Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat
Usage: python3 analysis/supernova_tide_test.py <Pantheon+SH0ES.dat> [--json out.json]
"""
import json
import sys

import numpy as np

C = 299792.458
f = sys.argv[1]
hdr = open(f).readline().split()
raw = np.genfromtxt(f, names=hdr, dtype=None, encoding=None, skip_header=1)
keep = (raw["zHD"] > 0.023) & (raw["IS_CALIBRATOR"] == 0)  # skip very nearby (own motions) and calibrators
d = raw[keep]
# the same supernova can appear from several surveys: combine them into one weighted point
by = {}
for r in d:
    by.setdefault(r["CID"], []).append(r)
z, m, e, ra, dec = [], [], [], [], []
for cid, rs in by.items():
    w = np.array([1 / r["m_b_corr_err_DIAG"] ** 2 for r in rs])
    z.append(np.mean([r["zHD"] for r in rs])); m.append(np.sum(w * np.array([r["m_b_corr"] for r in rs])) / w.sum())
    e.append(1 / np.sqrt(w.sum())); ra.append(rs[0]["RA"]); dec.append(rs[0]["DEC"])
z, m, e = np.array(z), np.array(m), np.array(e)
ra, dec = np.radians(ra), np.radians(dec)
nvec = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
w = 1 / e**2
print(f"supernovae used: {len(z)} (z = {z.min():.3f} to {z.max():.2f})")

# distance-modulus table mu[om, ol, sn] for each supernova's redshift (H0 absorbed into the offset)
OM = np.round(np.arange(0.0, 0.61, 0.025), 4)
OL = np.round(np.arange(-2.5, 2.51, 0.025), 4)
zg = np.linspace(0, z.max() * 1.01, 1500)


def dl(om, ol):
    ok = 1 - om - ol
    E2 = om * (1 + zg) ** 3 + ok * (1 + zg) ** 2 + ol
    if np.any(E2 <= 0):
        return None  # this universe re-collapses or never reaches these redshifts
    inv = 1 / np.sqrt(E2)
    chi = np.concatenate([[0], np.cumsum((inv[1:] + inv[:-1]) / 2 * np.diff(zg))])
    if ok > 1e-6:
        s = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok)
    elif ok < -1e-6:
        s = np.sin(np.sqrt(-ok) * chi) / np.sqrt(-ok)
    else:
        s = chi
    return (1 + zg) * s


MU = np.full((len(OM), len(OL), len(z)), np.nan)
for i, om in enumerate(OM):
    for j, ol in enumerate(OL):
        D = dl(om, ol)
        if D is not None:
            MU[i, j] = 5 * np.log10(np.interp(z, zg, D) * C)


def chi2_of(mu):
    """chi-square with the brightness/H0 offset fitted analytically; nan if the model can't reach some supernovae."""
    r = m - mu
    off = np.sum(w * r, axis=-1, keepdims=True) / w.sum()
    return np.sum(w * (r - off) ** 2, axis=-1)


def mu_dir(i, ol_sn):
    """distance moduli for Omega_m index i and a per-supernova Omega_L (linear interpolation in the table)."""
    x = (ol_sn - OL[0]) / (OL[1] - OL[0])
    j0 = np.clip(np.floor(x).astype(int), 0, len(OL) - 2)
    t = x - j0
    idx = np.arange(len(z))
    a, b = MU[i, j0, idx], MU[i, j0 + 1, idx]
    return a * (1 - t) + b * t


res = {}
# matter only
best = min((chi2_of(MU[i, int(np.argmin(np.abs(OL)))]), OM[i]) for i in range(len(OM)))
res["matter only"] = dict(chi2=float(best[0]), params=dict(Omega_m=float(best[1])), k=2)
# flat dark energy
cands = []
for i, om in enumerate(OM):
    j = int(np.argmin(np.abs(OL - (1 - om))))
    c = chi2_of(MU[i, j])
    if np.isfinite(c):
        cands.append((c, om))
best = min(cands)
res["dark energy (same everywhere)"] = dict(chi2=float(best[0]), params=dict(Omega_m=float(best[1])), k=2)

# axes: a Fibonacci half-sky (an axis and its opposite are the same for a tide)
NA = 600
k_ = np.arange(NA) + 0.5
zc = k_ / NA
phi = np.pi * (1 + 5 ** 0.5) * k_
axes = np.stack([np.sqrt(1 - zc**2) * np.cos(phi), np.sqrt(1 - zc**2) * np.sin(phi), zc], 1)
EPS = np.round(np.arange(-2.0, 2.01, 0.05), 3)


def tide_fit(ol0_of_om):
    best = (np.inf, None)
    P2 = (3 * (nvec @ axes.T) ** 2 - 1) / 2  # [sn, axis]
    for i, om in enumerate(OM):
        ol0 = ol0_of_om(om)
        for eps in EPS:
            ol_sn = ol0 + eps * P2.T  # [axis, sn]
            mu = np.stack([mu_dir(i, ol_sn[a]) for a in range(NA)])
            c = chi2_of(mu)
            c = np.where(np.isfinite(c), c, np.inf)
            a = int(np.argmin(c))
            if c[a] < best[0]:
                v = axes[a]
                l_, b_ = np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2]))
                best = (float(c[a]), dict(Omega_m=float(om), Omega_L0=float(ol0), eps=float(eps), axis_ra=float(l_), axis_dec=float(b_)))
    return best


# A big universe close to a well is pulled harder on the near side than the far side, so the stretch is lopsided:
# add a one-sided (dipole) term d * P1 on the same axis. The axis now needs the full sky (the sign matters).
FULL = np.concatenate([axes, -axes])
DIP = np.round(np.arange(-1.5, 1.51, 0.1), 3)
QUAD = np.round(np.arange(-2.0, 2.01, 0.1), 3)


def lopsided_fit(om, ol0):
    i = int(np.argmin(np.abs(OM - om)))
    cosang = nvec @ FULL.T                       # [sn, axis]
    P1, P2 = cosang.T, ((3 * cosang**2 - 1) / 2).T
    best = (np.inf, None)
    for dd in DIP:
        for qq in QUAD:
            mu = mu_dir_many(i, ol0 + dd * P1 + qq * P2)
            c = chi2_of(mu); c = np.where(np.isfinite(c), c, np.inf)
            a = int(np.argmin(c))
            if c[a] < best[0]:
                v = FULL[a]
                best = (float(c[a]), dict(Omega_m=float(OM[i]), Omega_L0=float(ol0), dipole=float(dd), quadrupole=float(qq),
                                          axis_ra=float(np.degrees(np.arctan2(v[1], v[0])) % 360), axis_dec=float(np.degrees(np.arcsin(v[2])))))
    return best


def mu_dir_many(i, ol):  # ol: [axis, sn]
    x = (ol - OL[0]) / (OL[1] - OL[0])
    j0 = np.clip(np.floor(x).astype(int), 0, len(OL) - 2)
    t = x - j0
    idx = np.broadcast_to(np.arange(len(z)), ol.shape)
    return MU[i, j0, idx] * (1 - t) + MU[i, j0 + 1, idx] * t


LOPSIDED = "--lopsided" in sys.argv
c, p = tide_fit(lambda om: 0.0)
res["tide only (no dark energy)"] = dict(chi2=c, params=p, k=5)
c, p = tide_fit(lambda om: 1 - om)
res["dark energy + tide"] = dict(chi2=c, params=p, k=5)
if LOPSIDED:
    for om in (0.3,):
        c, p = lopsided_fit(om, 0.0)
        res[f"lopsided well alone, matter {om}"] = dict(chi2=c, params=p, k=6)
        c, p = lopsided_fit(om, 1 - om)
        res[f"dark energy + lopsided pull, matter {om}"] = dict(chi2=c, params=p, k=6)
    # reference at the same matter content
    i = int(np.argmin(np.abs(OM - 0.3))); j = int(np.argmin(np.abs(OL - 0.7)))
    res["dark energy, matter 0.3"] = dict(chi2=float(chi2_of(MU[i, j])), params=dict(Omega_m=0.3), k=2)

n = len(z)
print(f"\n{'model':32s} {'chi-square':>11s} {'settings':>9s} {'BIC':>9s}")
base = res["dark energy (same everywhere)"]["chi2"] + res["dark energy (same everywhere)"]["k"] * np.log(n)
for name, r in res.items():
    bic = r["chi2"] + r["k"] * np.log(n)
    r["bic"] = float(bic); r["delta_bic_vs_dark_energy"] = float(bic - base)
    print(f"{name:32s} {r['chi2']:11.1f} {r['k']:9d} {bic:9.1f}   (vs dark energy: {bic - base:+.1f})   {r['params']}")
print("\nBIC = chi-square + (number of settings) x ln(number of supernovae); lower is better, a gap above ~10 is decisive.")
if "--json" in sys.argv:
    json.dump(dict(n=n, results=res), open(sys.argv[sys.argv.index("--json") + 1], "w"))
