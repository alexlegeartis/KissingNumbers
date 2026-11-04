import json,numpy as np
octads=json.load(open('research/mixed_tier17_20261104/octads18.json'))['octads'];masks=[sum(1<<j for j in s) for s in octads];rm1=sorted(set([0,65535]+masks));rm2=[x for x in range(65536) if all((x&o).bit_count()%2==0 for o in masks)];supports=[x for x in rm2 if x.bit_count()==6];assert len(rm1)==32 and len(rm2)==2048 and len(supports)==448
unseen=set(rm2);cosets=[]
while unseen:
 r=min(unseen);f={r^x for x in rm1};assert f<=unseen;unseen-=f;cosets.append(r)
assert len(cosets)==64
V=[];keys=[]
for s in supports:
 for parity in (0,1):
  # RM(1,4) acts transitively on each parity class on a bent 6-set.
  orbit={h&s for h in rm1};assert len(orbit)==32 and all(x.bit_count()%2==0 for x in orbit)
  sign=0 if parity==0 else s&-s;v=[(-2 if sign>>i&1 else 2) if s>>i&1 else 0 for i in range(16)];V.append(v);keys.append((0,s,sign))
for q in range(16):
 for w in cosets:
  v=[(-1 if w>>i&1 else 1)*(3 if i==q else 1) for i in range(16)];V.append(v);keys.append((1,q,w))
assert len(V)==1920
np.array(V,np.int8).tofile('research/mixed_tier17_20261104/mixed17_orbits.bin');np.array(rm1,np.uint16).tofile('research/mixed_tier17_20261104/mixed17_group.bin');json.dump(keys,open('research/mixed_tier17_20261104/mixed17_keys.json','w'));print('orbits',len(V),'supports',len(supports),'rm2 cosets',len(cosets))
np.array(cosets,np.uint16).tofile('research/mixed_tier17_20261104/mixed17_flat_reps.bin')
