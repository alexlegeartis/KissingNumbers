import numpy as np
from pathlib import Path
Path('research/mixed_tier17_order8_20261105/cache').mkdir(parents=True,exist_ok=True)
p='research/mixed_tier17_eval_subgroup_20261105/cache/';r='research/mixed_tier17_eval_subgroup_20261105/seed/';out='research/mixed_tier17_order8_20261105/cache/'
V=np.fromfile(p+'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile(p+'mixed17_eval4_group.bin',np.uint16);F=np.fromfile(p+'mixed17_eval4_flat_reps.bin',np.uint16);H0=np.fromfile(r+'mixed17_rank4_group.bin',np.uint16);T=np.fromfile(r+'mixed17_seed_tier.bin',np.int8).reshape(-1,16).astype(np.int16);flats=np.fromfile(r+'mixed17_seed_flatminus.bin',np.uint16)
S=lambda H:np.array([[1-2*(int(h)>>i&1) for i in range(16)] for h in H],np.int16)
neg_t={tuple(y) for t in T for y in t*S(H0)};neg_f={int(w)^int(h) for w in flats for h in H0};neg=[]
for i,w in enumerate(F):
 orb={int(w)^int(h) for h in H}
 if orb<=neg_f:neg.append(i);neg_f-=orb
for i,v in enumerate(V):
 orb={tuple(y) for y in v*S(H)}
 if orb<=neg_t:neg.append(128+i);neg_t-=orb
assert not neg_t and not neg_f and len(neg)==44
np.savetxt(out+'negative_seed.txt',neg,fmt='%d');print('reconstructed',len(neg),'negative order-16 seed orbits')
