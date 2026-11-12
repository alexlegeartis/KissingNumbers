import os,sys,json,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from scipy.optimize import milp,LinearConstraint,Bounds
from scipy.sparse import coo_matrix
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent
p=ROOT/'verifications/improved/dim41-46-p48-heads'
e=np.load(p/'data/e8_geometry.npz');u=e['u'];norm=e['norm'];ids=e['ids'];tiny=np.load(p/'data/population_table.npz')['tiny']
weights=np.zeros(len(u),dtype=np.int64)
weights[norm==0]=1092000
weights[norm==8]=116235
weights[norm==16]=9180
weights[norm==24]=495
weights[norm==32]=tiny[ids[norm==32]]
roots=u[norm==8];u=np.r_[u,3*roots];weights=np.r_[weights,np.ones(len(roots),dtype=np.int64)]
assert weights.sum()==52416000
# Simple roots: standard chain A; D chain ending in branch; E8 basis.
def normals(kind,k):
 if kind=='A':
  V=np.zeros((k,8),dtype=int)
  for i in range(k):V[i,i]=2;V[i,i+1]=-2
  return V
 if kind=='D':
  V=normals('A',k-1);v=np.zeros(8,dtype=int);v[k-2:k]=2
  return np.r_[V,[v]]
 if kind=='E':
  V=np.array([[1,-1,-1,-1,-1,-1,-1,1],[2,2,0,0,0,0,0,0],[-2,2,0,0,0,0,0,0],[0,-2,2,0,0,0,0,0],[0,0,-2,2,0,0,0,0],[0,0,0,-2,2,0,0,0],[0,0,0,0,-2,2,0,0],[0,0,0,0,0,-2,2,0]])
  return V[:k]
 raise ValueError(kind)

def run(kind,k):
 V=normals(kind,k);G=3*V@V.T/4;GI=np.linalg.inv(G)
 assert np.all((u@V.T)%4==0)
 profiles=u@V.T//4
 cells,inv=np.unique(profiles,axis=0,return_inverse=True);w=np.bincount(inv,weights=weights).astype(np.int64)
 norms=np.einsum('ij,jk,ik->i',cells,GI,cells)
 keep=(w>0)&(norms<6-1e-10);cells=cells[keep];w=w[keep];norms=norms[keep]
 dots=cells@GI@cells.T
 compat=4*(3-dots)**2<=(6-norms[:,None])*(6-norms)+1e-9
 compat|=dots>=3
 np.fill_diagonal(compat,False)
 # global antipodal parity of counts allows fixing an arbitrary orientation only for exact symmetry; do not assume.
 print(kind,k,'cells',len(cells),'edges',compat.sum()//2,'maxweight',max(w),'central',w[np.all(cells==0,axis=1)].tolist(),flush=True)
 # weighted greedy multistart; exact geometry checked later
 neigh=[set(np.flatnonzero(row)) for row in compat]
 best=[];bestw=0
 rng=np.random.default_rng(20261113)
 starts=list(range(len(cells)))
 for rep in range(5):
  rng.shuffle(starts)
  for s in starts:
   cand=neigh[s].copy();cl=[s]
   while cand:
    if rep==0:q=max(cand,key=lambda a:(w[a],len(neigh[a]&cand)))
    elif rep==1:q=max(cand,key=lambda a:w[a]*(1+len(neigh[a]&cand))**.2)
    else:q=max(cand,key=lambda a:w[a]*rng.uniform(.6,1.4))
    cl.append(q);cand&=neigh[q]
   ww=int(sum(w[cl]))
   if ww>bestw:bestw=ww;best=cl;print('best',kind,k,bestw,[cells[i].tolist() for i in best],flush=True)
  if len(cells)>3000:break
 out=dict(kind=kind,k=k,dimension=48-k,bound=bestw,normals=V.tolist(),gram=G.astype(int).tolist(),offsets=cells[best].tolist(),weights=w[best].tolist(),cells=len(cells),status='numeric geometry needs exact replay')
 (OUT/f'profile_{kind}{k}.json').write_text(json.dumps(out,indent=2))
 np.savez_compressed(OUT/f'profile_{kind}{k}.npz',cells=cells,weights=w,compat=compat,normals=V)
 return out
if __name__=='__main__':
 for spec in sys.argv[1:] or ['A1','A2','A3','D4','D5','E6','E7','E8']:
  run(spec[0],int(spec[1:]))
