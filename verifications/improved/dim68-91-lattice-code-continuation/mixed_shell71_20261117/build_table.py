import numpy as np,json
from pathlib import Path
P=Path(__file__).resolve().parent;D=np.load(P/'leech_structure.npz');A=D['a'];B=D['b'];F=D['gram'];assert np.array_equal(A@B,2*np.eye(24,dtype=int));assert np.array_equal(A@F@A.T,2*F);assert np.array_equal(B@F@B.T,2*F)
keys,cnts=np.unique(np.fromfile(P/'min4keys.bin',dtype=np.uint32),return_counts=True);assert len(keys)==98280 and np.all(cnts==2)
bkeys=np.array([sum((int(x)&1)<<i for i,x in enumerate(row))for row in B],dtype=np.uint32)
groups=np.zeros(len(keys),dtype=np.uint32)
for i in range(24):groups^=((keys>>i)&1)*bkeys[i]
ug,gc=np.unique(groups,return_counts=True);assert len(ug)==4095 and np.all(gc==24) and 0 not in ug
N=np.zeros(1<<24,dtype=np.uint16)
for g in ug:
 a=keys[groups==g];c=(a[:,None]^a[None,:]^g).ravel();np.add.at(N,c,1)
assert N.sum()==4095*576
N*=4
# kernelB=imageA generated overF2
span={0}
for row in A:
 key=sum((int(x)&1)<<i for i,x in enumerate(row));span|={s^key for s in list(span)}
assert len(span)==4096
for k in span:
 if k:N[k]+=96
assert N[0]==0 and N.sum()==9828000
N.tofile(P/'fiber8_table.bin');json.dump({'type2_cosets':len(keys),'frame_cosets':len(ug),'table_sum':int(N.sum()),'table_hist':{str(x):int(c)for x,c in zip(*np.unique(N,return_counts=True))}},open(P/'table_report.json','w'),indent=2)
print(open(P/'table_report.json').read())
