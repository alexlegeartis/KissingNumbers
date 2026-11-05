"""Independently enumerate the exact 33D one-coordinate support neighborhoods."""
import collections,itertools,json,importlib.util
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[2]/'verifications/improved/dim32-35-signed-layers/source/classic'
spec=importlib.util.spec_from_file_location('classic',root/'verify.py');classic=importlib.util.module_from_spec(spec);spec.loader.exec_module(classic)
p=root/'certificate33.json';assert classic.verify(p)['verified_points']==370368
F=json.loads(p.read_text())['full_even_supports'];assert len(F)==1848 and len(set(F))==1848
base=set(F);all_candidates=set();domain=(1<<33)-1
for t in F:
 on=t
 while on:
  bit=on&-on;on-=bit;off=domain^t
  while off:
   new=off&-off;off-=new;s=t^bit^new
   if s not in base:all_candidates.add(s)
assert len(all_candidates)==369600
F=np.asarray(F,dtype=np.uint64);S=np.array(sorted(all_candidates),dtype=np.uint64);hist=np.zeros(4,np.int64);groups=collections.defaultdict(list)
lut=np.array([i.bit_count() for i in range(256)],dtype=np.uint8)
def popcount(a):
 if hasattr(np,'bitwise_count'):return np.bitwise_count(a)
 return lut[a.view(np.uint8).reshape(a.shape+(8,))].sum(axis=-1,dtype=np.uint8)
for st in range(0,len(S),2048):
 part=S[st:st+2048];hit=popcount(part[:,None]&F[None,:])>=5;counts=hit.sum(axis=1);hist+=np.bincount(np.minimum(counts,3),minlength=4)
 for i in np.flatnonzero(counts==1):groups[int(np.flatnonzero(hit[i])[0])].append(int(part[i]))
assert hist.tolist()==[0,131,662,368807] and len(groups)==94
counts=collections.Counter()
for i,L in groups.items():
 union=int(F[i])
 for s in L:union|=s
 counts[len(L),union.bit_count()]+=1
expected={(1,9):69,(2,9):15,(2,10):3,(3,9):2,(4,9):5}
assert counts==expected and sum(v for (a,b),v in counts.items() if b==9)==91
remaining=[]
for i,L in groups.items():
 union=int(F[i]);common=union
 for s in L:union|=s;common&=s
 if union.bit_count()==10:
  assert common.bit_count()==6;remaining.append(i)
assert len(remaining)==3
print('PASS: 369600 candidate supports; 91 groups lie in nine coordinates and the remaining 3 share six coordinates')
