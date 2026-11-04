import os
import numpy as np,itertools,json,time,random
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
V=np.fromfile('research/mixed_tier17_20261104/mixed17_rank4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16);H=np.fromfile('research/mixed_tier17_20261104/mixed17_rank4_group.bin',np.uint16);T=np.fromfile('research/mixed_tier17_20261104/mixed17_rank4_same.bin',np.uint8).reshape(3840,128);G=np.fromfile('research/mixed_tier17_20261104/mixed17_rank4_flat_graph.bin',np.uint8).reshape(128,128);signs=np.array([[(-1 if int(h)>>k&1 else 1) for k in range(16)] for h in H],np.int16);baseline=[0,12,22,27,39,83,98,111,116,120,41,53];rng=random.Random(171141);seen=set();best=44;start=time.monotonic()
for trial in range(int(os.environ.get("MIXED17_TRIALS", "1500"))):
 flat=set(baseline);remove=rng.randint(1,9)
 for i in rng.sample(sorted(flat),remove):flat.remove(i)
 order=list(range(128));rng.shuffle(order);goal=rng.choice([9,10,11,12,13,14])
 for a in order:
  if len(flat)>=goal:break
  if a not in flat and all(not G[a,b] for b in flat):flat.add(a)
 if len(flat)<8:continue
 S=tuple(sorted(flat))
 if S in seen:continue
 seen.add(S);I=np.flatnonzero(np.all(T[:,S]==0,axis=1));X=V[I];n=len(I)
 if n+len(S)<=best:continue
 adj=np.zeros((n,n),bool)
 for sg in signs:adj|=((X*sg)@X.T)>8
 np.fill_diagonal(adj,False);rr,cc=np.nonzero(np.triu(adj,1));rows=np.repeat(np.arange(len(rr)),2);cols=np.column_stack((rr,cc)).ravel();A=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(rr),n)).tocsr();r=milp(-np.ones(n),integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,1),options={'time_limit':float(os.environ.get('MIXED17_MILP_SECONDS', '3'))})
 val=None if r.fun is None else int(round(-r.fun));upper=None if r.mip_dual_bound is None else -r.mip_dual_bound;total=0 if val is None else val+len(S)
 if total>best:
  best=total;chosen=I[np.flatnonzero(r.x>.5)];np.savez('research/mixed_tier17_20261104/mixed17_rank4_joint_candidate.npz',flat=S,tier=chosen);print('POTENTIAL_ONE_HEMISPHERE',json.dumps({'trial':trial,'flat':S,'tier':val,'total_orbits':total,'seconds':round(time.monotonic()-start,1)}),flush=True);break
 if trial%50==0:print(json.dumps({'trial':trial,'flat':len(S),'candidates':n,'incumbent':total,'upper':upper,'distinct':len(seen),'seconds':round(time.monotonic()-start,1)}),flush=True)
print('done best',best,'distinct',len(seen),flush=True)
