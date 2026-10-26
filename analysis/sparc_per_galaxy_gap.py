import json, glob, os, sys, numpy as np
"""Fit every SPARC galaxy on its own with a stacked hidden sheet, letting the gap between the sheets, the hidden mass,
its spread and the stars' mass-to-light ratio all be free. Then ask whether one common gap fits most galaxies.
Usage: python3 analysis/sparc_per_galaxy_gap.py <sparc_dir> [out.json]   (data as in sparc_tests.py)"""
D = sys.argv[1]; G = 4.30091e-6
meta = {}
for line in open(D + "/SPARC_Lelli2016c.mrt"):
    t = line.split()
    if len(t) >= 18 and t[1].isdigit():
        try: meta[t[0]] = dict(inc=float(t[5]), L=float(t[7]) * 1e9, rd=float(t[11]), mhi=float(t[13]) * 1e9, q=int(t[17]))
        except ValueError: pass
F = np.logspace(-0.5, 2.0, 26); K = np.logspace(-0.3, 1.5, 19); U = np.array([0.3, 0.4, 0.5, 0.65, 0.8])  # star mass-to-light, per galaxy
H = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 13])
rows = []
for f in sorted(glob.glob(D + "/rotmod/*_rotmod.dat")):
    name = os.path.basename(f).replace("_rotmod.dat", ""); m = meta.get(name)
    if not m or m["q"] >= 3 or m["inc"] < 30: continue
    a = np.loadtxt(f, comments="#"); r, vo, ev, vg, vd, vb = a.T[:6]
    ok = (r > 0) & (vo > 0); r, vo, ev, vg, vd, vb = r[ok], vo[ok], ev[ok], vg[ok], vd[ok], vb[ok]
    if len(r) < 5: continue
    sig = np.sqrt(np.maximum(2 * ev / (vo * np.log(10)), 0.02) ** 2 + 0.05 ** 2); lo = np.log10(vo ** 2 / r); rd = max(m["rd"], 0.1)
    prof = []
    for h in H:
        rho = np.sqrt(r ** 2 + h ** 2); best = np.inf
        for u in U:
            gb = (vg * np.abs(vg) + u * vd * np.abs(vd) + (u + 0.2) * vb * np.abs(vb)) / r
            if np.any(gb <= 0): gb = np.maximum(gb, 1e-3)
            mb = u * m["L"] + 1.33 * m["mhi"]
            for k in K:
                shape = G * mb / (rho + k * rd) ** 2 * (r / rho)
                c2 = (((np.log10(gb[None, :] + F[:, None] * shape[None, :]) - lo) / sig) ** 2).sum(1).min()
                best = min(best, c2)
        prof.append(best)
    prof = np.array(prof); i = int(np.argmin(prof))
    scale = max(1.0, prof[i] / max(1, len(r) - 4))  # widen the allowed range when the fit is worse than the error bars say
    okh = H[prof - prof[i] <= scale]
    rows.append(dict(name=name, rd=rd, h=float(H[i]), lo=float(okh.min()), hi=float(okh.max()), dchi0=float(prof[0] - prof[i]), redchi=float(prof[i] / max(1, len(r) - 4))))
R = rows; con = [x for x in R if x["hi"] - x["lo"] <= 3 and x["hi"] < 13]
h = np.array([x["h"] for x in con]); rd = np.array([x["rd"] for x in con]); mm = h > 0
print(f"galaxies {len(R)}, gap pinned down in {len(con)}; pinned gaps median {np.median(h):.2f} kpc, 16-84%: {np.percentile(h,16):.2f}..{np.percentile(h,84):.2f}")
print(f"gap vs disc size slope {np.polyfit(np.log10(rd[mm]), np.log10(h[mm]), 1)[0]:.2f} (0 = universal, 1 = scales), correlation {np.corrcoef(np.log10(rd[mm]), np.log10(h[mm]))[0,1]:.2f}")
for hc in [0, 1, 2, 3, 4]:
    print(f"  one common gap {hc} kpc consistent with {100*np.mean([x['lo'] <= hc <= x['hi'] for x in R]):.0f}% of galaxies")
print("galaxies clearly preferring some gap over none:", round(float(np.mean([x['dchi0'] > 4 for x in R])), 2))
print("median fit quality (reduced chi-square):", round(float(np.median([x['redchi'] for x in R])), 2))
if len(sys.argv) > 2: json.dump(R, open(sys.argv[2], "w"))
