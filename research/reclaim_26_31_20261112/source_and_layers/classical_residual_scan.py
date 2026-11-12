"""Probe exact sqrt(3)/2 Leech lifts against current PackingStar coordinates."""
import argparse, itertools, json, math, time
from pathlib import Path
import numpy as np
from source_paths import PACKINGSTAR_DEFAULT, packingstar_file, fresh_leech_shell

p=argparse.ArgumentParser();p.add_argument('--dim',type=int,required=True);p.add_argument('--packingstar',type=Path,default=PACKINGSTAR_DEFAULT);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
d=args.dim;k=d-24;t0=time.time()
root=args.output;root.mkdir(parents=True,exist_ok=True)
source=packingstar_file(args.packingstar,d)
f=np.load(source,allow_pickle=False)
name='points' if 'points' in f else 'rational_points'
den=int(f['coordinate_denominator' if 'coordinate_denominator' in f else 'rational_denominator'])
x=f[name].astype(np.int16)
assert np.all(np.sum(x.astype(np.int64)**2,axis=1)==32*den*den)
L=fresh_leech_shell()
core=x[np.all(x[:,24:]==0,axis=1),:24]
assert np.all(core%den==0)
core_set=set(map(bytes,(core//den).astype(np.int16)))
assert core_set <= set(map(bytes,L))
q=L[[bytes(v) not in core_set for v in L]]
tails=[]
for a,b in itertools.combinations(range(k),2):
 for sa,sb in itertools.product([-2,2],repeat=2):
  tt=np.zeros(k,dtype=np.int16);tt[a]=sa;tt[b]=sb;tails.append(tt)
tails=np.array(tails)
allowed=np.ones((len(q),len(tails)),dtype=bool)
tailvals,inv,cnts=np.unique(x[:,24:],axis=0,return_inverse=True,return_counts=True)
order=sorted(np.flatnonzero(np.any(tailvals!=0,axis=1)),key=lambda j:-cnts[j])
print('start',d,'missing',len(q),'tailvalues',len(tails),'groups',len(order),flush=True)
for z,g in enumerate(order):
 ids=np.flatnonzero(allowed.any(axis=1))
 if len(ids)==0:break
 a=x[inv==g,:24].astype(np.float32); b=tailvals[g]
 # All dot products here are sums of small exactly represented integers.
 # Cauchy bounds every partial absolute sum by sqrt(32*32)*den <=64.
 bound=np.array([math.isqrt((4*int(16*den-tt@b)**2)//3) for tt in tails],dtype=np.int16)
 for begin in range(0,len(ids),512):
  jj=ids[begin:begin+512]
  maxima=(q[jj].astype(np.float32)@a.T).max(axis=1)
  allowed[jj]&=maxima[:,None]<=bound[None,:]
 if z%10==0 or z+1==len(order):print('group',z+1,'of',len(order),'surviving',int(allowed.sum()),'heads',int(allowed.any(1).sum()),'sec',time.time()-t0,flush=True)

ids,jj=np.nonzero(allowed); cq=q[ids];ct=tails[jj]
if 'quadratic_A' in f:
 A=f['quadratic_A'];B=f['quadratic_B'];ds=f['quadratic_denominators'];rs=f['quadratic_radicands']
 assert np.all(A[:,:24]==0) and np.all(B[:,24:]==0) and np.all(ds==1) and np.all(rs==3)
 keep=np.all(3*cq.astype(np.int32)@B[:,:24].T+2*ct.astype(np.int32)@A[:,24:].T<=32,axis=1)
 # Equality with an already present point fails its self inner product test.
 cq,ct=cq[keep],ct[keep]

out=root/f'classical_residual_{d}.npz';np.savez_compressed(out,leech_heads=cq,tails=ct)
report={'dim':d,'missing_Leech_vectors':len(q),'tail_directions':len(tails),'initial_candidates':len(q)*len(tails),'candidate_count':len(cq),'seconds':time.time()-t0,'file':str(out),'scope':'Fixed-source classical sqrt(3)/2 Leech lifts only; no global exclusion'}
(root/f'classical_residual_{d}.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
