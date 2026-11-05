from pathlib import Path
Path('research/mixed_tier17_eval_subgroup_20261105/cache').mkdir(parents=True,exist_ok=True)
import json,numpy as np
o=json.load(open('research/mixed_tier17_20261104/octads18.json'))['octads'];octads=[sum(1<<i for i in a) for a in o];full=sorted(set([0,65535]+octads));H=[h for h in full if (h&1)==0];assert len(H)==16
rm2=[w for w in range(65536) if all((w&o).bit_count()%2==0 for o in octads)];supp=[w for w in rm2 if w.bit_count()==6];unseen=set(rm2);flat=[]
while unseen:
 r=min(unseen);f={r^h for h in H};assert f<=unseen;unseen-=f;flat.append(r)
assert len(flat)==128
V=[];keys=[]
for s in supp:
 proj={h&s for h in H};assert len(proj)==16;all_sign={w for w in range(65536) if w&~s==0}
 while all_sign:
  w=min(all_sign);orb={w^x for x in proj};assert orb<=all_sign;all_sign-=orb;V.append([(-2 if w>>i&1 else 2) if s>>i&1 else 0 for i in range(16)]);keys.append((0,s,w))
for q in range(16):
 for w in flat:V.append([(-1 if w>>i&1 else 1)*(3 if i==q else 1) for i in range(16)]);keys.append((1,q,w))
assert len(V)==3840
np.asarray(V,np.int8).tofile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_orbits.bin');np.asarray(H,np.uint16).tofile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin');np.asarray(flat,np.uint16).tofile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_flat_reps.bin');json.dump(keys,open('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_keys.json','w'));print('tier orbits',len(V),'flat orbits',len(flat))

np.asarray([0,15,20,27,38,41,51,60,71,72,81,94,141,142,145,146,297,298,309,310,709,710,729,730,1053,1054,1065,1066,1141,1142,1153,1154,1489,1490,1509,1510,1593,1594,1613,1614,1697,1698,1725,1726],np.int32).tofile('research/mixed_tier17_eval_subgroup_20261105/cache/baseline_seed.bin')
