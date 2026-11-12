"""Probe exact sqrt(3)/2 Leech lifts against current PackingStar coordinates."""
import argparse, itertools, json, math, time
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--dim',type=int,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--leech',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
d=args.dim;k=d-24;t0=time.time()
root=Path(__file__).parent
f=np.load(args.source,allow_pickle=False)
name='points' if 'points' in f else 'rational_points'
den=int(f['coordinate_denominator' if 'coordinate_denominator' in f else 'rational_denominator'])
x=f[name].astype(np.int16)
if x.shape[1] != d or not np.all(np.sum(x.astype(np.int64)**2,axis=1)==32*den*den):
 raise ValueError('Source dimension or normalization')
L=np.load(args.leech,allow_pickle=False).astype(np.int16)
if L.shape != (196560,24) or not np.all(np.sum(L.astype(np.int64)**2,axis=1)==32):
 raise ValueError('Leech shell normalization')
core=x[np.all(x[:,24:]==0,axis=1),:24]
if np.any(core%den):raise ValueError('Core is not an integral Leech subset')
core_set=set(map(bytes,(core//den).astype(np.int16)))
if not core_set <= set(map(bytes,L)):
 raise ValueError('Classical completion lemma requires an equatorial Leech subset')
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
 thresholds=[int(16*den-tt@b) for tt in tails]
 if min(thresholds)<0:raise ValueError('Unexpected negative radical threshold')
 bound=np.array([math.isqrt((4*v*v)//3) for v in thresholds],dtype=np.int16)
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

out=args.output;np.savez_compressed(out,leech_heads=cq,tails=ct)
print(json.dumps({'dim':d,'candidate_count':len(cq),'seconds':time.time()-t0,'file':str(out)}),flush=True)
