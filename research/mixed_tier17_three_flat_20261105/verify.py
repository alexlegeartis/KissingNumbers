"""Exact three-flat-coset obstruction. Run from the repository root."""
import json,itertools
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
octads=json.loads((root/'research/mixed_tier17_20261104/octads18.json').read_text())['octads'];masks=[sum(1<<j for j in s) for s in octads];H=sorted(set([0,65535]+masks));assert len(H)==32 and all((x^y) in H for x in H for y in H)
rm2=[x for x in range(65536) if all((x&o).bit_count()%2==0 for o in masks)];supports=[x for x in rm2 if x.bit_count()==6];assert len(rm2)==2048 and len(supports)==448
unseen=set(rm2);F=[]
while unseen:
 r=min(unseen);coset={r^x for x in H};assert coset<=unseen;unseen-=coset;F.append(r)
assert len(F)==64
V=[]
for s in supports:
 orbit={h&s for h in H};assert len(orbit)==32 and all(x.bit_count()%2==0 for x in orbit)
 for parity in (0,1):
  sign=0 if parity==0 else s&-s;V.append([(-2 if sign>>i&1 else 2) if s>>i&1 else 0 for i in range(16)])
for j in range(16):
 for word in F:V.append([(-1 if word>>i&1 else 1)*(3 if i==j else 1) for i in range(16)])
V=np.asarray(V,dtype=np.int16);H=np.asarray(H,dtype=np.int64);F=np.asarray(F,dtype=np.int64);assert V.shape==(1920,16)
signs=np.array([[(-1 if int(h)>>i&1 else 1) for i in range(16)] for h in H],np.int16);words=np.array([[(-1 if (int(f)^int(h))>>i&1 else 1) for i in range(16)] for f in F for h in H],np.int16)
same=(V@words.T).reshape(1920,64,32).max(axis=2)>8
flat=np.zeros((64,64),bool)
for i in range(64):
 for j in range(i):flat[i,j]=flat[j,i]=min((int(F[i])^int(F[j])^int(h)).bit_count() for h in H)<6
powers=1<<np.arange(16,dtype=np.int64);support=((V!=0)*powers).sum(1);negative=((V<0)*powers).sum(1);cs=np.min(negative[:,None]^(support[:,None]&H[None,:]),axis=1);keys={(tuple(abs(v)),int(t)):i for i,(v,t) in enumerate(zip(V,cs))};assert len(keys)==1920
canon=lambda x:min(int(x)^int(h) for h in H);fk={canon(v):i for i,v in enumerate(F)};assert len(fk)==64

def tiermap(W):
 supp=((W!=0)*powers).sum(1);neg=((W<0)*powers).sum(1);c=np.min(neg[:,None]^(supp[:,None]&H[None,:]),axis=1)
 return np.array([keys[(tuple(abs(v)),int(t))] for v,t in zip(W,c)])
def maskperm(f,p):return sum(((int(f)>>p[x])&1)<<x for x in range(16))
def check_action(fm,tm):
 assert len(set(fm))==64 and len(set(tm))==1920
 assert np.array_equal(flat,flat[np.ix_(fm,fm)]) and np.array_equal(same,same[np.ix_(tm,fm)])
 # The transformation is a signed coordinate permutation preserving H, so it
 # preserves every tier-tier conflict exactly, not just the seed covers.

shears=[]
for i in range(4):
 for j in range(4):
  if i==j:continue
  p=np.array([x^(((x>>j)&1)<<i) for x in range(16)]);assert {maskperm(h,p) for h in H}==set(H)
  fm=np.array([fk[canon(maskperm(f,p))] for f in F]);tm=tiermap(V[:,p]);check_action(fm,tm);shears.append(fm)
trans=[]
for f in F:
 sg=np.array([-1 if int(f)>>i&1 else 1 for i in range(16)]);fm=np.array([fk[canon(int(x)^int(f))] for x in F]);tm=tiermap(V*sg);check_action(fm,tm);trans.append(fm)
nb=[i for i in range(1,64) if not flat[0,i]];expected={(0,)+c for c in itertools.combinations(nb,2) if not any(flat[a,b] for a,b in itertools.combinations(c,2))};assert len(expected)==168
cert=json.loads((here/'fractional_cover.json').read_text());S=tuple(cert['flat']);cliques=cert['cliques'];weights=cert['weights'];den=cert['denominator'];assert type(den) is int and den>0 and len(cliques)==len(weights) and all(type(w) is int and w>0 for w in weights)
allowed=set(np.flatnonzero(np.all(~same[:,S],axis=1)).tolist());cover={v:0 for v in allowed}
for block,weight in zip(cliques,weights):
 assert len(block)==len(set(block)) and set(block)<=allowed
 for a,b in itertools.combinations(block,2):assert int(np.max(signs@(V[a]*V[b])))>8
 for v in block:cover[v]+=weight
assert all(x>=den for x in cover.values()) and sum(weights)<20*den
orbit={S};queue=[S]
for T in queue:
 for fm in shears:
  U=[int(fm[x]) for x in T]
  for z in U:
   W=tuple(sorted(int(trans[z][x]) for x in U))
   if W not in orbit:orbit.add(W);queue.append(W)
assert orbit==expected
print('verified 168 normalized three-flat sets by exact symmetry')
print('verified rational clique cover weight',sum(weights),'/',den,'< 20')
print('PASS: the tier has at most 19 orbits, so three flats plus tier has at most 22 orbits.')
