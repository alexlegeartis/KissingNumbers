import numpy as np,json
R='research/mixed_tier17_20261104/';C=np.fromfile(R+'mixed17_same.bin',np.uint8).reshape(1920,64);Splus=[0,12,22,27,39,41];Sminus=[1,13,23,26,38,40];K=json.load(open(R+'mixed17_keys.json'));KR=json.load(open(R+'mixed17_rank4_keys.json'));full=set(map(int,np.fromfile(R+'mixed17_group.bin',np.uint16)));F=np.fromfile(R+'mixed17_flat_reps.bin',np.uint16);FR=np.fromfile(R+'mixed17_rank4_flat_reps.bin',np.uint16);T=np.fromfile(R+'mixed17_rank4_orbits.bin',np.int8).reshape(-1,16);parents=np.flatnonzero(np.all(C[:896,Splus]==0,axis=1));assert len(parents)==16 and not np.any(C[np.ix_(parents,Sminus)])
short={(K[int(i)][1],K[int(i)][2].bit_count()%2) for i in parents};tier=[i for i,(kind,s,w) in enumerate(KR) if kind==0 and (s,w.bit_count()%2) in short];assert len(tier)==32
def flats(S):
 return [i for i,r in enumerate(FR) if any((int(r)^int(F[j])) in full for j in S)]
p=flats(Splus);m=flats(Sminus);assert len(p)==len(m)==12 and not set(p)&set(m)
np.asarray(T[tier],np.int8).tofile(R+'mixed17_seed_tier.bin');np.asarray(FR[p],np.uint16).tofile(R+'mixed17_seed_flatplus.bin');np.asarray(FR[m],np.uint16).tofile(R+'mixed17_seed_flatminus.bin');print('seed tier orbits',len(tier),'flat halves',len(p),len(m),'baseline',4322+2*(32+12)*16)
