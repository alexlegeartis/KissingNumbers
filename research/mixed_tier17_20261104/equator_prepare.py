import json,itertools,numpy as np
R='research/mixed_tier17_20261104/';octads=json.load(open(R+'octads18.json'))['octads'];O=[sum(1<<i for i in a) for a in octads];H=list(map(int,np.fromfile(R+'mixed17_rank4_group.bin',np.uint16)));V=[];meta=[]
def rep_orbits(mask,words):
 projection={h&mask for h in H};pending=set(words)
 while pending:
  w=min(pending);orbit={w^h for h in projection};assert orbit<=pending;pending-=orbit;yield w,len(orbit)
for a,b in itertools.combinations(range(16),2):
 s=(1<<a)|(1<<b);words=[0,1<<a,1<<b,s]
 for w,n in rep_orbits(s,words):V.append([(-4 if w>>i&1 else 4) if s>>i&1 else 0 for i in range(16)]);meta.append((s,w,n,0))
for s in O:
 positions=[i for i in range(16) if s>>i&1];words=[sum(((x>>j)&1)<<i for j,i in enumerate(positions)) for x in range(256) if x.bit_count()%2]
 for w,n in rep_orbits(s,words):V.append([(-2 if w>>i&1 else 2) if s>>i&1 else 0 for i in range(16)]);meta.append((s,w,n,1))
assert sum(a[2] for a in meta)==4320
X=[];info=[]
for positions in itertools.combinations(range(16),6):
 s=sum(1<<j for j in positions);words=[sum(((x>>j)&1)<<i for j,i in enumerate(positions)) for x in range(64)]
 for w,n in rep_orbits(s,words):X.append([(-2 if w>>i&1 else 2) if s>>i&1 else 0 for i in range(16)]);info.append((s,w,n))
print('equator orbits',len(V),'candidate orbits',len(X),'candidate sizes',dict(zip(*np.unique(np.array(info)[:,2],return_counts=True))),flush=True)
np.asarray(V,np.int8).tofile(R+'mixed17_eq_orbits.bin');np.asarray(meta,np.int32).tofile(R+'mixed17_eq_meta.bin');np.asarray(X,np.int8).tofile(R+'mixed17_six_all.bin');np.asarray(info,np.int32).tofile(R+'mixed17_six_meta.bin')
