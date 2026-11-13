import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,itertools,sys,json
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
O=Path(__file__).parent;R=O.parent/'kn_latest'
sys.path.insert(0,str(R/'research/reclaim_26_31_20261113/algebraic'))
from optimize_ls15_orientation import extract
p=np.load(O.parent/'design31_20261113/design.npz')['points'];rng=np.random.default_rng(92830)
for dim,count in [(30,248084)]:
 x=np.load(R/f'verifications/improved/dim30-leech-sloane-section/data/coordinates{dim}_{count}.npz')['points'].astype(np.int16)//2
 E=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,dim)]; xp=x[:,E];ids=np.flatnonzero(np.sum(xp.astype(int)**2,axis=1)>8);old=xp[ids].astype(np.float32);keys={r.tobytes() for r in x}
 best=0;trials=0
 for omitted in [()]:
  base=np.load(O/'native30module.npz')['points']
  base=base[np.sum(base[:,8:].astype(int)**2,axis=1)==16]
  y=base.copy()
  for trial in range(1):
   for _ in range(0 if trial==0 else 4):
    mode=rng.integers(4)
    if mode==0:
     i,j=rng.choice(8,2,False);y[:,[i,j]]=y[:,[j,i]]
    elif mode==1:
     i,j=8+rng.choice(dim-24,2,False);y[:,[i,j]]=y[:,[j,i]]
    elif mode==2:
     i,j=rng.choice(8,2,False);y[:,[i,j]]*=-1
    else:
     v=rng.choice([-1,1],8);v[-1]=np.prod(v[:-1]);d=y[:,:8]@v;assert not np.any(d%4);y[:,:8]-=d[:,None]//4*v
   full=np.zeros((len(y),dim),np.int16);full[:,E]=y
   new=np.array([r for r in full if r.tobytes() not in keys],dtype=np.int16)
   if not len(new):continue
   g=csr_matrix(new[:,E].astype(np.float32)@old.T>16);m=maximum_bipartite_matching(g,perm_type='column');gain=len(new)-sum(m>=0);trials+=1
   if gain>best:
    best=int(gain);ai,ri=extract(g,m);keep=np.ones(len(x),bool);keep[ids[ri]]=False;f=np.r_[x[keep],new[ai]]
    for s in range(0,len(ai),64):
     dot=new[ai[s:s+64]].astype(np.float64)@x[keep].T.astype(np.float64);assert dot.max(initial=-32)<=16
    dot=new[ai].astype(int)@new[ai].T.astype(int);np.fill_diagonal(dot,-32);assert dot.max(initial=-32)<=16
    np.savez_compressed(O/f'candidate{dim}_{len(f)}.npz',points=f);np.savez_compressed(O/f'delta{dim}_{len(f)}.npz',added=new[ai],removed=x[ids[ri]])
    print('BEST',dim,omitted,trial,best,'total',len(f),flush=True)
  print('section',dim,omitted,'heads',len(base),'best_gain',best,flush=True)
 (O/f'report{dim}.json').write_text(json.dumps(dict(dimension=dim,trials=trials,best_gain=best)))
