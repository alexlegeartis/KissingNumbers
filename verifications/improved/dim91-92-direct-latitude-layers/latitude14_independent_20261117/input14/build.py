from pathlib import Path
import numpy as np,json,hashlib
P=Path(__file__).resolve().parent
S=P.parent/'hexagon_sparse15_phase7812_20261117'
A=np.load(P.parent/'pless48_fiber9_cert_20261117/min15.npy')
G=np.load(P.parent/'pless48_fiber9_cert_20261117/generator.npy'); np.save(P/'generator.npy',G)
K=np.load(S/'phase_code.npy'); np.save(P/'source_phase7812.npy',K)
n,j,v=max((int(np.sum(K[:,j]==v)),j,v) for j in range(15) for v in range(3))
H=np.delete(K[K[:,j]==v],j,axis=1);np.save(P/'phase14.npy',H)
hist=np.zeros(15,dtype=np.int64)
for i in range(len(H)-1): hist+=np.bincount(np.count_nonzero(H[i+1:]!=H[i],axis=1),minlength=15)
assert hist[:5].sum()==0 and hist[5]>0
assert np.array_equal(G[:,:24],np.eye(24,dtype=np.uint8))
rows=[]
for i in range(3):
 for t in (1,2):
  p=np.zeros(3,np.uint8);p[i]=t
  full=A[np.all(A[:,:3]==p,axis=1)]
  assert len(full)==31680 and np.all(np.count_nonzero(full[:,3:],axis=1)==14)
  assert np.array_equal((full[:,:24].astype(np.int16)@G.astype(np.int16))%3,full)
  assert len(np.unique(full,axis=0))==len(full)
  fn='outer_'+''.join(map(str,p))+'.npy';np.save(P/fn,full[:,3:])
  rows.append(dict(prefix=p.tolist(),file=fn,count=len(full)))
report=dict(phase_count=n,shortening_coordinate=j,shortening_value=v,phase_distance_histogram=hist.tolist(),outer_layers=rows,points_per_layer=31680*n,addition91=2*31680*n,addition92=4*31680*n,equator=21943640934,total91_without_poles=21943640934+2*31680*n,total92_without_poles=21943640934+4*31680*n,total91=21943640936+2*31680*n,total92=21943640940+4*31680*n)
(P/'certificate.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
