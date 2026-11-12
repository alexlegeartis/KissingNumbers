import os
os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['OMP_NUM_THREADS']='2'
import json,time,sys,itertools
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).parent;P=ROOT/'verifications/improved/dim49-60-p48-pointwise-refill'
a=np.load(P/'heads49.npz');G=a['gram'].astype(np.int64)
orig=np.r_[a['base'],a['extra_points']].astype(np.int64)
raw=np.load(P/'discovery/low_conflict_pool.npz')['points'].astype(np.int64)
raw=np.unique(np.r_[raw,-raw,orig,-orig],axis=0)
S=np.load(OUT/'p48_signed_class.npz')['points']
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 20261113)

def clique(adj,target,cand=None):
 if target==0:return []
 if cand is None:cand=(1<<len(adj))-1
 if cand.bit_count()<target:return None
 while cand:
  bit=cand&-cand;v=bit.bit_length()-1;cand^=bit
  if cand.bit_count()+1<target:return None
  z=clique(adj,target-1,cand&adj[v])
  if z is not None:return [v]+z
 return None

history=[]
for rnd in range(20):
 start=time.time();keys=set(map(tuple,S));C=np.array([x for x in raw if tuple(x) not in keys]);CG=C@G
 assert 48*int(abs(CG).max())*int(abs(S).max())<2**53
 blocks=[];groups=defaultdict(list)
 for i in range(0,len(C),512):
  dots=CG[i:i+512].astype(float)@S.astype(float).T
  for q,row in enumerate(dots):
   b=tuple(map(int,np.flatnonzero(row>2)));blocks.append(b)
   if len(b)<=6:groups[b].append(i+q)
 print('round',rnd,'size',len(S),'hist',dict(Counter(map(len,blocks))),flush=True)
 active=np.ones(len(S),dtype=bool);chosen=[];deleted=[]
 # For each direct blocker set up to size6, all candidate subsets are retrieved.
 bs=list(groups.keys());rng.shuffle(bs);bs.sort(key=lambda b:len(b))
 for b in bs:
  k=len(b)
  if k>6 or any(not active[j] for j in b):continue
  js=[]
  for n in range(k+1):
   for sub in itertools.combinations(b,n):js.extend(groups.get(sub,()))
  if len(js)<=k:continue
  if chosen:
   js=[j for j in js if max(CG[j]@C[chosen].T)<=2]
  if len(js)<=k:continue
  points=C[js];dots=points@G@points.T
  adj=[]
  for i,row in enumerate(dots):
   flags=(row<=2);flags[i]=False
   adj.append(int.from_bytes(np.packbits(flags,bitorder='little').tobytes(),'little'))
  sel=clique(adj,k+1)
  if sel is None:continue
  picked=[js[j] for j in sel]
  active[list(b)]=False;chosen.extend(picked);deleted.extend(b)
  print('trade',k,'for',len(sel),'gain so far',len(chosen)-len(deleted),flush=True)
 if not chosen:break
 old=len(S);S=np.r_[S[active],C[chosen]]
 history.append(dict(round=rnd,before=old,after=len(S),removed=len(deleted),added=len(chosen),elapsed=time.time()-start))
 np.savez_compressed(OUT/'p48_signed_class_lns.npz',points=S,gram=G)
 (OUT/'p48_lns.json').write_text(json.dumps(dict(count=len(S),baseline=15598,gain=len(S)-15598,history=history),indent=2))
 print('END',len(S),'gain',len(S)-15598,'seconds',time.time()-start,flush=True)
# Verify geometry; integer products via exact binary64 with proven range.
assert len(set(map(tuple,S)))==len(S)
assert np.all(np.einsum('ij,ij->i',S@G,S)==6)
assert 48*int(abs(S@G).max())*int(abs(S).max())<2**53
for i in range(0,len(S),256):
 dots=(S[i:i+256]@G).astype(float)@S.astype(float).T
 for j in range(len(dots)):dots[j,i+j]=-6
 assert dots.max()<=2
keys=set(map(tuple,S));print('PASS',len(S),'asymmetric',sum(tuple(-x) not in keys for x in S),flush=True)
