"""Exact bounded-integer algebraic proposals from the small40D affine shell."""
from pathlib import Path
import argparse,json,numpy as np,time
from math import isqrt
ap=argparse.ArgumentParser();ap.add_argument('--inner',type=int,default=12);args=ap.parse_args();inner=args.inner;norm=60-2*inner
limits=[isqrt(n*norm)//2 for n in (54,36,30)]
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1];SOURCE=ROOT/'verifications/improved/dim40-affine-shells/source'
Q=np.array(json.loads((SOURCE/'gram.json').read_text()),np.int64)
z=np.load(SOURCE/'auxiliary.npz');B=z['norm36'].astype(np.int64);C=z['norm30'].astype(np.int64)
A=np.load(SOURCE/'core.npz')['coefficients']
if 40*40*int(abs(Q).max())*int(abs(C).max())**2>=2**62:raise ValueError('Source Gram integer bound')
start=time.monotonic();G=(C@Q)@C.T
if not np.all(np.diag(G)==30):raise ValueError('Incorrect source auxiliary norm')
u,c=np.unique(G,return_counts=True);print('C Gram',dict(zip(map(int,u),map(int,c))),flush=True)
rows,cols=np.where(np.triu(G==inner,1));D=C[rows]-C[cols];D=np.unique(np.r_[D,-D],axis=0)
print('difference norm',norm,'proposals',len(D),flush=True)
rng=np.random.default_rng(20261104);sample=A[np.sort(rng.choice(len(A),4096,replace=False))].astype(np.int64)
# Sparse norm36 candidates, all retained point decisions exact integers in
# binary32 under the proved absolute accumulation bound below2^24.
def dots(X,Y):
 if 40*int(abs(X).max(initial=0))*int(abs(Q).max(initial=0))>=2**62:raise ValueError('Coefficient/Gram integer bound')
 left=X@Q
 bound=int(abs(left).max(initial=0))*int(np.abs(Y.astype(np.int64)).sum(axis=1).max(initial=0))
 if bound>=2**24:raise ValueError(('binary32 exactness bound',bound))
 return left.astype(np.float32)@Y.astype(np.float32).T
keep=[];hist={};lows=[]
for at in range(0,len(D),256):
 x=D[at:at+256]
 mx=dots(x,sample).max(axis=1)
 good=x[mx<=limits[0]]
 if len(good):
  badB=(dots(good,B)>limits[1]).sum(axis=1);badC=(dots(good,C)>limits[2]).sum(axis=1)
  for v,b,cost in zip(good,badB,badC):
   key=(int(b),int(cost));hist[key]=hist.get(key,0)+1
  mask=badB==0;keep.extend(good[mask].tolist());lows.extend(map(int,badC[mask]))
print('sampleA compatible and B,C blocker histogram',sorted(hist.items())[:30],'B-free',len(keep),'min C blockers',min(lows,default=None),flush=True)
if keep:
 X=np.array(keep,np.int64);valid=np.ones(len(X),bool);maxima=np.full(len(X),-10**10,dtype=np.int64)
 for at in range(0,len(A),8192):
  ids=np.flatnonzero(valid)
  if not len(ids):break
  values=dots(X[ids],A[at:at+8192]).max(axis=1)
  maxima[ids]=np.maximum(maxima[ids],values.astype(np.int64));valid[ids[values>limits[0]]]=False
 print('all-core and B compatible',int(valid.sum()),flush=True)
 np.savez_compressed(OUT/f'difference40_q{inner}_candidates.npz',candidates=X[valid],core_maxima=maxima[valid],gram=Q)
 report={'inner':inner,'candidate_norm':norm,'limits':limits,'generated_differences':len(D),'sample_core_points':len(sample),'sample_and_B_survivors':len(keep),'full_core_and_B_survivors':int(valid.sum()),'C_blocker_counts':[int(x) for x in (dots(X[valid],C)>limits[2]).sum(axis=1)],'seconds':time.monotonic()-start,'new_sota':False,'scope':'All-core/B compatibility exact under explicit binary32 accumulation bounds. Candidate/C and candidate/candidate assembly still required.'}
else:report={'inner':inner,'candidate_norm':norm,'generated_differences':len(D),'sample_and_B_survivors':0,'seconds':time.monotonic()-start,'new_sota':False}
(OUT/f'difference40_q{inner}_probe.json').write_text(json.dumps(report,indent=2)+'\n')
print(report,flush=True)
