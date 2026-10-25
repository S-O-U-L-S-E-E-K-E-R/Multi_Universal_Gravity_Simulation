import numpy as np
G=6.674e-11; Msun=1.989e30; kpc=3.086e19; Gly=9.461e24
# galaxy: Hernquist visible mass
Mvis=6e10*Msun; a=3*kpc
def v_vis(r): return np.sqrt(G*Mvis*r/(r+a)**2)
rng=np.random.default_rng(1)
Muni=1.5e53  # kg, ordinary matter in observable universe ~1.5e53
# 30 neighbour universes at random directions, 60-200 Gly
N=30; dirs=rng.normal(size=(N,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
D=rng.uniform(60,200,N)*Gly; pos=dirs*D[:,None]; M=np.full(N,Muni)
def a_ext(x):
    d=pos[None,:,:]-x[:,None,:]; r=np.linalg.norm(d,axis=2)
    return (G*M[None,:,None]*d/r[:,:,None]**3).sum(1)
for rk in [5,10,20,30,50]:
    r=rk*kpc; th=np.linspace(0,2*np.pi,360,endpoint=False)
    x=np.stack([r*np.cos(th),r*np.sin(th),0*th],1)
    tid=a_ext(x)-a_ext(np.zeros((1,3)))  # relative to galaxy centre
    ar=-(tid*x).sum(1)/r   # inward positive
    need=(220e3**2/r)-v_vis(r)**2/r
    print(f"r={rk}kpc vis={v_vis(r)/1e3:.0f}km/s  needed inward extra={need:.2e}  tidal mean={ar.mean():.2e} min={ar.min():.2e} max={ar.max():.2e}")
# Gauss flux through sphere r=30kpc
r=30*kpc; u=rng.normal(size=(20000,3)); u/=np.linalg.norm(u,axis=1)[:,None]
tid=a_ext(u*r)-a_ext(np.zeros((1,3)))
print("mean inward tidal over sphere:", (-(tid*u).sum(1)).mean(), " (Gauss => 0)")
# mass needed at 100 Gly to give 1.5e-10 inward at 30kpc on the best axis (perp. axis gives -GMr/D^3 i.e. inward)
need=1.5e-10; Dn=100*Gly
Mneed=need*Dn**3/(G*r); print("M needed at 100Gly:", Mneed, "kg =", Mneed/Muni, "observable universes")
