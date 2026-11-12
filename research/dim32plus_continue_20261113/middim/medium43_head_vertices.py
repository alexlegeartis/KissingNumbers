"""New head templates from the fixed large-family polytope.

Only the large families and pure heads are fixed here. Medium and tiny
points may be removed by the later exact exchange computation.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,sys,math
import numpy as np
from scipy.spatial import HalfspaceIntersection
from medium43_pool import W,P,HERE,cutoff
sys.path.insert(0,str(W));import verify_families as vf

def run(family=59):
    g,_,_,_=vf.geometry();f=np.load(W/'configuration_43.npz');r=12
    C=f['classes'];R=f['r3'];D=12000000000000;den=int(f['head_denominator']);Z=f['head_numerators'].astype(float)/den
    pure=f['pure_directions'].astype(float)
    dc=((g['digits'][C]-g['digits'][int(C[family])])%3)@g['powers']
    beta=r+R-18+g['mu'][dc]//4
    ns=R+np.sum(Z*Z,axis=1);major=R!=12
    planes=np.c_[beta[major]/np.sqrt(r*ns[major]),2*Z[major]/np.sqrt(ns[major,None]),-np.ones(sum(major))]
    planes=np.r_[planes,np.c_[np.zeros(len(pure)),2*pure/np.linalg.norm(pure,axis=1)[:,None],-np.ones(len(pure))],[[1,0,0,0,-1],[-1,0,0,0,0]]]
    V=HalfspaceIntersection(planes,np.array([.01,0,0,0]),qhull_options='QJ1e-10').intersections
    headlist=[];Zi=f['head_numerators'].astype(object)*(D//den);NN=R.astype(object)*D*D+np.sum(Zi*Zi,axis=1)
    for v in V:
        if np.linalg.norm(v)<=1+1e-8 or v[0]<=1e-7:continue
        h=np.sqrt(r)*v[1:]/v[0]
        if np.max(np.abs(h))*D>2**62:continue
        z=tuple(map(int,np.rint(h*D)));n=r*D*D+sum(a*a for a in z)
        if any(int(beta[j])>cutoff(z,n,Zi[j],NN[j],D) for j in np.flatnonzero(major)):continue
        if any((w:=2*sum(a*int(b) for a,b in zip(z,p)))>0 and w*w>n*sum(int(b)*int(b) for b in p) for p in f['pure_directions']):continue
        if z not in headlist:headlist.append(z)
    result=dict(family=family,head_denominator=D,vertices_total=len(V),heads=headlist,status='Every head exactly checked against fixed major families and pure heads; medium/tiny deletions unresolved')
    (HERE/f'medium43_extraheads_f{family}.json').write_text(json.dumps(result,indent=2))
    print('VERTICES',len(V),'HEADS',len(headlist),'radii',[round(sum(a*a for a in z)/D**2,4) for z in headlist],flush=True)
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 59)
