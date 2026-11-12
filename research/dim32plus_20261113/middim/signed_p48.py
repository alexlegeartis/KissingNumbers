import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import sys,time,json
from pathlib import Path
import numpy as np
from collections import Counter,defaultdict
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;P=ROOT/'verifications/improved/dim49-60-p48-pointwise-refill'
a=np.load(P/'heads49.npz');G=a['gram'].astype(np.int64);base=np.r_[a['base'],a['extra_points']].astype(np.int64)
raw=np.load(P/'discovery/low_conflict_pool.npz')['points'].astype(np.int64)
# Canonicalize candidate lines, drop already selected lines.
def canon(x):
 x=x.copy();sgn=np.where(x[np.arange(len(x)),np.argmax(x!=0,axis=1)]<0,-1,1);return x*sgn[:,None]
base=canon(base);raw=canon(raw)
basekeys=set(map(tuple,base));raw=np.array([r for r in raw if tuple(r) not in basekeys],dtype=np.int64)
raw=np.unique(raw,axis=0)
S=np.r_[base,-base];C=np.r_[raw,-raw]
assert np.all(np.einsum('ij,ij->i',C@G,C)==6)
print('baseline',len(S),'candidates',len(C),flush=True)
# Exact binary64 integer products: absolute accumulation bound.
assert 48*int(abs(C@G).max())*int(abs(S).max())<2**53
sf=S.astype(np.float64);cg=(C@G).astype(np.float64)
blockers=[]
for i in range(0,len(C),512):
 dots=cg[i:i+512]@sf.T
 for row in dots:blockers.append(set(map(int,np.flatnonzero(row>2))))
print('hist',dict(Counter(map(len,blockers))),flush=True)
active=np.ones(len(S),dtype=bool);chosen=[];trades=[]
rng=np.random.default_rng(20261113)

def refresh():
 groups=defaultdict(list);free=[]
 if chosen:
  dots=cg@C[chosen].astype(float).T
  blocked=(dots>2).any(axis=1)
 else:blocked=np.zeros(len(C),bool)
 chosenset=set(chosen)
 for i,b in enumerate(blockers):
  if blocked[i] or i in chosenset:continue
  b.intersection_update(set(np.flatnonzero(active))) if False else None
  bb=tuple(sorted(j for j in b if active[j]))
  if len(bb)<=2:groups[bb].append(i)
 return groups

for rnd in range(100):
 groups=refresh();change=False
 free=groups.get((),[]);rng.shuffle(free)
 for i in free:
  if not chosen or np.all((C[chosen]@G@C[i])<=2):chosen.append(i);change=True;trades.append(dict(remove=[],add=[i]))
 # Each proposal is checked against all current additions; selected old rows covered by blockers.
 singles=[(b,js) for b,js in groups.items() if len(b)==1 and len(js)>=2];rng.shuffle(singles)
 for b,js in singles:
  if not active[b[0]]:continue
  js=[i for i in js if (not chosen or np.all(C[chosen]@G@C[i]<=2))]
  if len(js)<2:continue
  cc=C[js];dots=cc@G@cc.T
  found=False
  for ii in range(len(js)):
   jj=np.flatnonzero(dots[ii,ii+1:]<=2)
   if len(jj):
    sel=[js[ii],js[ii+1+int(jj[0])]];active[b[0]]=False;chosen.extend(sel);trades.append(dict(remove=list(b),add=sel));change=True;found=True;break
 if not change:break
 print('round',rnd,'selected',int(active.sum())+len(chosen),'gain',int(active.sum())+len(chosen)-len(S),'trades',len(trades),flush=True)
 final=np.r_[S[active],C[chosen]]
 np.savez_compressed(OUT/'p48_signed_class.npz',points=final,gram=G,removed=np.flatnonzero(~active),added=C[chosen])
 (OUT/'p48_signed_search.json').write_text(json.dumps(dict(baseline=len(S),count=len(final),gain=len(final)-len(S),rounds=rnd+1,trades=trades),indent=2))
# Full integer-coordinate independent validation in binary64 exact arithmetic.
final=np.r_[S[active],C[chosen]]
assert len(set(map(tuple,final)))==len(final)
assert np.all(np.einsum('ij,ij->i',final@G,final)==6)
f=final.astype(float);fg=(final@G).astype(float)
assert 48*int(abs(final@G).max())*int(abs(final).max())<2**53
maxdot=-6
for i in range(0,len(final),256):
 dots=fg[i:i+256]@f.T
 for j in range(len(dots)):dots[j,i+j]=-6
 maxdot=max(maxdot,int(dots.max()))
 assert dots.max()<=2
print('EXACT PAIR PASS',len(final),'gain',len(final)-len(S),'maxdot',maxdot,flush=True)
