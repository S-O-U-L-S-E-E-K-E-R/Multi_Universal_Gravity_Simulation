"""Our universe on the slope of a neighbour's well: stretching by direction and the mark it leaves on the CMB.
For a neighbour as heavy as our observable universe at several distances, prints the change in expansion rate on the
near side, far side and sides, the leftover CMB pattern (Sachs-Wolfe estimate, order of magnitude) and, since both
scale with mass, the heaviest neighbour the CMB allows (pattern < 1e-5) and the mass needed to replace dark energy.
Usage: python3 analysis/neighbour_limits.py"""
import numpy as np
G=6.674e-11; c=2.998e8; Gly=9.461e24; H0=70e3/3.086e22; Muni=1.5e53; R=46.5*Gly
def sphere(n=2000):
    k=np.arange(n)+0.5; z=1-2*k/n; p=np.pi*(1+5**.5)*k; s=np.sqrt(1-z*z); return np.stack([s*np.cos(p),s*np.sin(p),z],1)
U=sphere()
def effects(M, Dr):  # neighbour mass (x our universe), distance in our radii, along +x
    P=np.array([Dr*R,0,0]); x=U*R
    d=P-x; acc=G*M*Muni*(d/np.linalg.norm(d,axis=1)[:,None]**3 - P/np.linalg.norm(P)**3)
    ar=(acc*U).sum(1)                       # outward relative pull at our horizon, each direction
    dH=ar/(H0**2*R)                         # change in expansion rate, as a fraction (over ~one Hubble time)
    phi=-G*M*Muni/np.linalg.norm(d,axis=1); g0=G*M*Muni*P/np.linalg.norm(P)**3
    tid=phi-(-G*M*Muni/np.linalg.norm(P))-(-(g0*x).sum(1))   # potential left after removing the uniform fall
    dT=(tid.max()-tid.min())/2/(3*c*c)       # Sachs-Wolfe size of the leftover pattern on the CMB
    return dH, dT
for Dr in [1.2,2,10,100,1000]:
    dH,dT=effects(1,Dr)
    # mass allowed by the CMB (pattern < 1e-5) and mass needed to replace dark energy (outward pull ~0.7 H0^2 R on the near side)
    Mcmb=1e-5/dT; Mde=0.7/dH.max()
    print(f"D={Dr:6g} R: neighbour as heavy as our universe -> near-side expansion +{100*dH.max():.2g}%, far {100*dH[np.argmin(U[:,0])]:.2g}%, sides {100*dH[np.argmin(abs(U[:,0]))]:.2g}%, CMB pattern {dT:.1e} | allowed mass < {Mcmb:.2g}x, needed for dark energy {Mde:.2g}x")
