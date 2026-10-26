import numpy as np
rng = np.random.default_rng(42)
def vec(l, b):
    l, b = np.radians(l), np.radians(b)
    return np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)])
def lb(v):
    v = v/np.linalg.norm(v); return (np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2])))
def axis_sep(a, b):  # angle between two axes (sign ignored), 0..90
    return np.degrees(np.arccos(min(1, abs(a @ b))))
anoms = {
  "Radio+quasar excess dipole (kinematic part removed)": (217, 20),   # Secrest et al. 2022
  "Bulk flow, CosmicFlows-4":                           (297, -4),   # Watkins et al. 2023
  "Hubble-rate anisotropy, X-ray clusters":             (280, -15),  # Migkas et al. 2021
  "CMB hemispherical power asymmetry":                  (221, -22),  # Planck 2020
  "CMB 'axis of evil' (quadrupole-octopole)":           (240, 63),   # de Oliveira-Costa 2004 / later ILC
  "Dark flow (cluster kSZ)":                            (290, 30),   # Kashlinsky et al.
}
ref = {"CMB dipole (our motion)": (264.0, 48.3), "Quasar dipole, raw (Secrest 2021)": (238.2, 28.8), "Radio dipole NVSS, raw": (196, 46)}
names = list(anoms); V = np.array([vec(*anoms[n]) for n in names])
def best_axis(V):
    T = sum(np.outer(v, v) for v in V); w, e = np.linalg.eigh(T); return e[:, -1]
def stat(V):  # mean angle of the axes from their best common axis
    a = best_axis(V); return np.mean([axis_sep(v, a) for v in V])
A = best_axis(V); s_obs = stat(V)
print("Best common axis (l, b) = (%.0f, %.0f) [or the opposite direction]" % lb(A))
for n in names: print(f"  {axis_sep(vec(*anoms[n]), A):5.1f} deg  from common axis : {n} {anoms[n]}")
for n, d in ref.items(): print(f"  {axis_sep(vec(*d), A):5.1f} deg  (reference) : {n} {d}")
print("Mean angle from common axis: %.1f deg" % s_obs)
def rand_axes(k):
    z = rng.uniform(-1, 1, k); p = rng.uniform(0, 2*np.pi, k); s = np.sqrt(1-z*z)
    return np.stack([s*np.cos(p), s*np.sin(p), z], 1)
sims = np.array([stat(rand_axes(len(V))) for _ in range(200000)])
p = np.mean(sims <= s_obs)
print("Random 6-axis sets: median mean-angle %.1f deg; fraction as tight or tighter than observed: p = %.4f" % (np.median(sims), p))
# leave-one-out robustness
for i, n in enumerate(names):
    Vi = np.delete(V, i, 0); si = stat(Vi); simsi = np.array([stat(rand_axes(len(Vi))) for _ in range(40000)])
    print(f"  without {n[:45]:45s}: mean angle {si:4.1f}, p = {np.mean(simsi <= si):.4f}")
# pairwise table
print("\nPairwise axis separations (deg):")
short = ["Excess", "BulkFl", "H0anis", "HPA", "AoE", "DarkFl"]
print("        " + " ".join(f"{s:>7s}" for s in short))
for i in range(len(V)): print(f"{short[i]:>7s} " + " ".join(f"{axis_sep(V[i], V[j]):7.0f}" for j in range(len(V))))

# The three motion anomalies (bulk flow, H0 anisotropy, dark flow) form a tight group. Is a group this tight among
# 6 random axes unusual? Search every triple in each random set, since we picked the triple after looking.
from itertools import combinations
idx = [names.index(n) for n in ["Bulk flow, CosmicFlows-4", "Hubble-rate anisotropy, X-ray clusters", "Dark flow (cluster kSZ)"]]
t_obs = stat(V[idx]); At = best_axis(V[idx])
print("\nMotion group (bulk flow, H0 anisotropy, dark flow): mean angle %.1f deg from axis (l,b)=(%.0f,%.0f)" % ((t_obs,) + lb(At if lb(At)[0] > 180 else -At)))
cnt = 0; N = 40000
for _ in range(N):
    R = rand_axes(6)
    if min(stat(R[list(c)]) for c in combinations(range(6), 3)) <= t_obs: cnt += 1
print("Chance that 6 random axes contain some triple this tight: p = %.3f" % (cnt / N))
print("Shapley Supercluster (l,b)~(306,30) is %.0f deg from that axis; Great Attractor (307,9) is %.0f deg" % (axis_sep(vec(306,30), At), axis_sep(vec(307,9), At)))
