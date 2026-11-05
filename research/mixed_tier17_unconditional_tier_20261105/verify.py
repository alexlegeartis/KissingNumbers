"""Exact structural checks for the global same-cap tier bound."""
from pathlib import Path
import json,itertools,subprocess,sys,numpy as np
root=Path(__file__).resolve().parents[2];octads=json.loads((root/'research/mixed_tier17_20261104/octads18.json').read_text())['octads'];masks=[sum(1<<j for j in s) for s in octads];H=sorted(set([0,65535]+masks));assert len(H)==32 and all((a^b) in H for a in H for b in H)
rm2=[x for x in range(65536) if all((x&o).bit_count()%2==0 for o in masks)];supports=[x for x in rm2 if x.bit_count()==6];assert len(rm2)==2048 and len(supports)==448
remaining=set(rm2);F=[]
while remaining:
 r=min(remaining);coset={r^h for h in H};assert coset<=remaining;remaining-=coset;F.append(r)
assert len(F)==64
short=[];labels=[]
for s in supports:
 orbit={h&s for h in H};assert len(orbit)==32 and all(x.bit_count()%2==0 for x in orbit)
 for parity in (0,1):
  sign=0 if parity==0 else s&-s;short.append([(-2 if sign>>i&1 else 2) if s>>i&1 else 0 for i in range(16)]);labels.append(s)
short=np.asarray(short,dtype=np.int16);labels=np.asarray(labels,dtype=np.uint16);signs=np.array([[(-1 if h>>i&1 else 1) for i in range(16)] for h in H],dtype=np.int16);assert short.shape==(896,16)
# All short-orbit pairs with >=3 shared support coordinates are incompatible.
conf=np.zeros((896,896),bool)
for sg in signs:conf|=(short*sg)@short.T>8
pop=np.array([i.bit_count() for i in range(65536)],dtype=np.uint8);need=pop[np.bitwise_and(labels[:,None],labels[None,:])]>=3;np.fill_diagonal(need,False);assert np.all(conf[need])
# Each dual coordinate gives a 64-orbit clique; any short orbit containing
# that coordinate conflicts with all of that dual-coordinate clique.
duals=[]
for q in range(16):
 dual=np.array([[(-1 if word>>i&1 else 1)*(3 if i==q else 1) for i in range(16)] for word in F],dtype=np.int16);cross=short[(labels>>q)&1==1];duconf=np.zeros((64,64),bool);crossconf=np.zeros((len(cross),64),bool)
 for sg in signs:
  duconf|=(dual*sg)@dual.T>8;crossconf|=(cross*sg)@dual.T>8
 np.fill_diagonal(duconf,True);assert np.all(duconf) and np.all(crossconf);duals.append(dual)
print('PASS: all support-intersection, dual-clique and support/dual conflicts hold exactly for the complete RM(1,4)-orbit family.')
print('The centered-support simplex argument in README gives <=16 short orbits, and then <=20 total short+dual orbits.')

# Complete the flat-count case analysis. All checks are exact integer checks.
V=np.r_[short,np.concatenate(duals)];words=np.array([[(-1 if (f^h)>>i&1 else 1) for i in range(16)] for f in F for h in H],np.int16);same=(V@words.T).reshape(1920,64,32).max(axis=2)>8
flat=np.zeros((64,64),np.int64)
for i in range(64):
 for j in range(i):flat[i,j]=flat[j,i]=int(min((F[i]^F[j]^h).bit_count() for h in H)<6)
assert np.all(flat.sum(axis=1)==35);I=np.eye(64,dtype=np.int64);assert not np.any((flat+5*I)@(flat-3*I)@(flat-35*I))
# Symmetric 35-regular graph with least eigenvalue >= -5: Hoffman gives
# alpha <= 64*5/(35+5)=8, so at most eight flat cosets can be chosen.
nb=[i for i in range(1,64) if not flat[0,i]]
cert=json.loads((root/'research/mixed_tier17_20261104/mixed17_five_flat_clique_certificate.json').read_text());expected={(0,)+c for c in itertools.combinations(nb,4) if not any(flat[a,b] for a,b in itertools.combinations(c,2))};assert len(expected)==280 and {tuple(x['flat']) for x in cert}==expected
for item in cert:
 S=tuple(item['flat']);eligible=set(np.flatnonzero(np.all(~same[:,S],axis=1)).tolist());cliques=item['cliques'];used=[v for c in cliques for v in c];assert len(cliques)==16 and len(used)==len(set(used)) and set(used)==eligible
 for c in cliques:
  for a,b in itertools.combinations(c,2):assert int(np.max(signs@(V[a]*V[b])))>8
# Every six-flat choice contains a five-flat subset and hence has <=16 tier orbits.
count=0
for c in itertools.combinations(nb,6):
 if any(flat[a,b] for a,b in itertools.combinations(c,2)):continue
 S=(0,)+c;eligible=np.flatnonzero(np.all(~same[:,S],axis=1));assert len(eligible)==16
 for a,b in itertools.combinations(eligible,2):assert int(np.max(signs@(V[a]*V[b])))>8
 count+=1
assert count==56
# The independent three/four-flat verifiers reconstruct their graphs and check
# their integer weighted/clique covers and exact symmetry actions.
for path in ('research/mixed_tier17_three_flat_20261105/verify.py','research/mixed_tier17_four_flat_20261105/verify.py'):
 subprocess.run([sys.executable,str(root/path)],cwd=root,check=True)
print('COMPLETE RESTRICTED FAMILY: k=0..2 use tier<=20; k=3 uses <=19; k=4 uses <=18; k=5,6 use <=16; k=7,8 use <=1. In every case k+tier<=22 orbits.')
