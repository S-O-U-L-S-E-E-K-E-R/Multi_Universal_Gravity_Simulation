"""Supernova test of the full picture: a lopsided pull toward a bigger well, a rotating universe, and our own motion.

Effective dark-energy term in each direction n (matter fixed at 30%):
  dark energy          : Omega_L(n) = M0
  lopsided well alone  : d * P1(n.a) + q * P2(n.a)                  (a pull from outside averages to zero)
  rotation alone       : W * sin^2(angle to spin axis)              (centrifugal push, zero along the axis)
  rotation + lopsided  : W * sin^2 + d * P1 + q * P2 on one axis     (W >= 0)

Two frames:
  A  standard: redshifts corrected to the CMB rest frame and for nearby flows (zHD)
  B  "we are moving, frame unknown": raw Sun-centred redshifts (zHEL) with our velocity fitted freely
     (direction and speed), because at universe scale our own motion can't be assumed.

Usage: python3 analysis/supernova_motion_rotation_test.py <Pantheon+SH0ES.dat> [--json out.json]
"""
import json
import sys

import numpy as np

C = 299792.458
OMEGA_M = 0.3
f = sys.argv[1]
hdr = open(f).readline().split()
raw = np.genfromtxt(f, names=hdr, dtype=None, encoding=None, skip_header=1)
d = raw[(raw["zHD"] > 0.023) & (raw["IS_CALIBRATOR"] == 0)]
by = {}
for r in d:
    by.setdefault(r["CID"], []).append(r)
zhd, zhel, m, e, ra, dec = [], [], [], [], [], []
for rs in by.values():
    w = np.array([1 / r["m_b_corr_err_DIAG"] ** 2 for r in rs])
    zhd.append(np.mean([r["zHD"] for r in rs])); zhel.append(np.mean([r["zHEL"] for r in rs]))
    m.append(np.sum(w * np.array([r["m_b_corr"] for r in rs])) / w.sum()); e.append(1 / np.sqrt(w.sum()))
    ra.append(rs[0]["RA"]); dec.append(rs[0]["DEC"])
zhd, zhel, m, e = map(np.array, (zhd, zhel, m, e))
ra, dec = np.radians(ra), np.radians(dec)
nvec = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
w = 1 / e**2
N = len(m)

# distance modulus table on (Omega_L, log z) at fixed matter
OL = np.round(np.arange(-3.0, 3.001, 0.025), 4)
LZ = np.linspace(np.log10(0.004), np.log10(2.5), 700)
zg = np.concatenate([[0], 10**LZ])
TAB = np.full((len(OL), len(LZ)), np.nan)
for j, ol in enumerate(OL):
    ok = 1 - OMEGA_M - ol
    E2 = OMEGA_M * (1 + zg) ** 3 + ok * (1 + zg) ** 2 + ol
    if np.any(E2 <= 0):
        continue
    inv = 1 / np.sqrt(E2)
    chi = np.concatenate([[0], np.cumsum((inv[1:] + inv[:-1]) / 2 * np.diff(zg))])
    s = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-6 else (np.sin(np.sqrt(-ok) * chi) / np.sqrt(-ok) if ok < -1e-6 else chi)
    dl = (1 + zg) * s
    with np.errstate(divide="ignore", invalid="ignore"):
        TAB[j] = 5 * np.log10(dl[1:] * C)


def mu(ol, z):
    """bilinear interpolation of the table; ol and z broadcast together."""
    x = (ol - OL[0]) / (OL[1] - OL[0]); j0 = np.clip(np.floor(x).astype(int), 0, len(OL) - 2); t = x - j0
    y = (np.log10(np.clip(z, 10**LZ[0], 10**LZ[-1])) - LZ[0]) / (LZ[1] - LZ[0]); k0 = np.clip(np.floor(y).astype(int), 0, len(LZ) - 2); u = y - k0
    return ((TAB[j0, k0] * (1 - u) + TAB[j0, k0 + 1] * u) * (1 - t) + (TAB[j0 + 1, k0] * (1 - u) + TAB[j0 + 1, k0 + 1] * u) * t)


def chi2(mu_):
    r = m - mu_
    off = np.sum(w * r, axis=-1, keepdims=True) / w.sum()
    c = np.sum(w * (r - off) ** 2, axis=-1)
    return np.where(np.isfinite(c), c, np.inf)


def sphere(nn, half=False):
    k = np.arange(nn) + 0.5
    zc = k / nn if half else 1 - 2 * k / nn
    ph = np.pi * (1 + 5**0.5) * k
    return np.stack([np.sqrt(1 - zc**2) * np.cos(ph), np.sqrt(1 - zc**2) * np.sin(ph), zc], 1)


AX = sphere(160)          # full sky (a one-sided pull has a direction); ~16 degrees apart
VDIR = sphere(300)
SPEEDS = np.arange(0, 1501, 50)


def radec(v):
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2])))


def model_ol(kind, p, cosang):
    if kind == "dark energy":
        return p["M0"] + 0 * cosang
    rot = p.get("W", 0.0) * (1 - cosang**2)
    lop = p.get("d", 0.0) * cosang + p.get("q", 0.0) * (3 * cosang**2 - 1) / 2
    return rot + lop


GRIDS = {
    "dark energy": dict(M0=np.round(np.arange(0.0, 1.41, 0.02), 3)),
    "lopsided well alone": dict(d=np.round(np.arange(-1.6, 1.61, 0.2), 2), q=np.round(np.arange(-2.4, 2.41, 0.3), 2)),
    "rotation alone": dict(W=np.round(np.arange(0.0, 2.51, 0.05), 3)),
    "rotation + lopsided well": dict(W=np.round(np.arange(0.0, 2.01, 0.2), 2), d=np.round(np.arange(-1.2, 1.21, 0.3), 2), q=np.round(np.arange(-1.6, 1.61, 0.4), 2)),
}
NSET = {"dark energy": 1, "lopsided well alone": 4, "rotation alone": 3, "rotation + lopsided well": 5}


def fit_model(kind, z):
    """best settings (and axis) for a given set of cosmological redshifts z."""
    grid = GRIDS[kind]; keys = list(grid)
    combos = np.array(np.meshgrid(*[grid[k] for k in keys], indexing="ij")).reshape(len(keys), -1).T
    cos_all = nvec @ AX.T  # [sn, axis]
    best = (np.inf, None, None)
    for combo in combos:
        p = dict(zip(keys, combo))
        if kind == "dark energy":
            c = chi2(mu(model_ol(kind, p, 0.0), z)); a = None
            if c < best[0]:
                best = (float(c), p, None)
            continue
        ol = model_ol(kind, p, cos_all.T)  # [axis, sn]
        c = chi2(mu(ol, z[None, :]))
        a = int(np.argmin(c))
        if c[a] < best[0]:
            best = (float(c[a]), p, AX[a])
    return best


def fit_velocity(kind, p, axis):
    """best observer velocity for given model settings, from Sun-centred redshifts."""
    cosang = nvec @ axis if axis is not None else 0.0
    ol = model_ol(kind, p, cosang)
    best = (np.inf, None)
    proj = nvec @ VDIR.T  # [sn, dir]
    for s in SPEEDS:
        zc = (1 + zhel[None, :]) / (1 - s * proj.T / C) - 1  # [dir, sn]
        c = chi2(mu(np.broadcast_to(ol, zc.shape), zc))
        a = int(np.argmin(c))
        if c[a] < best[0]:
            best = (float(c[a]), s * VDIR[a])
    return best


res = {"A_standard_frame": {}, "B_motion_free": {}}
for kind in GRIDS:
    c, p, ax = fit_model(kind, zhd)
    res["A_standard_frame"][kind] = dict(chi2=c, k=NSET[kind] + 1, params={k: float(v) for k, v in p.items()}, axis=radec(ax) if ax is not None else None)
    # frame B: alternate between our velocity and the model settings
    v = np.zeros(3)
    for it in range(2):
        cz = (1 + zhel) / (1 - (nvec @ v) / C) - 1
        c, p, ax = fit_model(kind, cz)
        cv, v = fit_velocity(kind, p, ax)
    res["B_motion_free"][kind] = dict(chi2=min(c, cv), k=NSET[kind] + 1 + 3, params={k: float(vv) for k, vv in p.items()}, axis=radec(ax) if ax is not None else None,
                                      our_velocity_kms=float(np.linalg.norm(v)), our_velocity_dir=radec(v / max(np.linalg.norm(v), 1e-9)))
    print(kind, "done", flush=True)

for frame, rr in res.items():
    base = rr["dark energy"]["chi2"] + rr["dark energy"]["k"] * np.log(N)
    print(f"\n{frame}  ({N} supernovae, matter {OMEGA_M})")
    for kind, r in rr.items():
        bic = r["chi2"] + r["k"] * np.log(N); r["delta_bic"] = float(bic - base)
        extra = f"  our velocity {r['our_velocity_kms']:.0f} km/s toward RA/Dec {r['our_velocity_dir'][0]:.0f},{r['our_velocity_dir'][1]:.0f}" if "our_velocity_kms" in r else ""
        print(f"  {kind:26s} chi2 {r['chi2']:7.1f}  vs dark energy BIC {bic - base:+6.1f}  {r['params']}  axis {r['axis']}{extra}")
print("\nFor reference, our motion measured from the CMB: 370 km/s toward RA/Dec 168,-7.")
if "--json" in sys.argv:
    json.dump(dict(n=N, results=res), open(sys.argv[sys.argv.index("--json") + 1], "w"))
