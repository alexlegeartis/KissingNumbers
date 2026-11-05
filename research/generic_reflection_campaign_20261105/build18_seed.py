from pathlib import Path
Path('research/generic_reflection_campaign_20261105/cache').mkdir(parents=True,exist_ok=True)
"""Build 18D bent-hexagon seed over Q(sqrt3) with metric diag(1^16,2,2)."""
import json,itertools,numpy as np
p='verifications/improved/dim18-bent-hexagon/data/';octads=json.load(open(p+'octads.json'))['octads'];fam=json.load(open(p+'families.json'))['families_at_angles_0_60_120'];tb=json.load(open(p+'tierB.json'))['slots_deg_to_words'];rows=[]
def signed(support,parity,scale):
 out=[]
 for bits in range(1<<len(support)):
  if bits.bit_count()%2!=parity:continue
  v=np.zeros(16,np.int64)
  for j,i in enumerate(support):v[i]=scale*(1-2*(bits>>j&1))
  out.append(v)
 return out
# Pairs (rational, radical) represent twice cos(k*pi/6).
cos=[(2,0),(0,1),(1,0),(0,0),(-1,0),(0,-1),(-2,0),(0,-1),(1*-1,0),(0,0),(1,0),(0,1)]
def add(v,L,k):
 ca,cb=cos[k%12];sa,sb=cos[(3-k)%12];a=np.r_[v,L*ca//2,L*sa//2];b=np.r_[np.zeros(16,dtype=np.int64),L*cb//2,L*sb//2];rows.append((a,b))
for i,j in itertools.combinations(range(16),2):
 for si,sj in itertools.product((-1,1),repeat=2):
  v=np.zeros(16,np.int64);v[i]=12*si;v[j]=12*sj;add(v,0,0)
for o in octads:
 for v in signed(o,1,6):add(v,0,0)
for pos in range(6):
 for support in fam[pos%3]:
  for v in signed(support,1,6):add(v,6,2*pos)
 add(np.zeros(16,np.int64),12,2*pos)
for deg,words in tb.items():
 for w in words:
  v=np.array([4*(1-2*((w>>(15-i))&1)) for i in range(16)],np.int64);add(v,4,int(deg)//30)
A=np.array([r[0] for r in rows]);B=np.array([r[1] for r in rows]);g=np.r_[np.ones(16,np.int64),[2,2]];assert A.shape==(8358,18);assert np.all(((A*A+3*B*B)*g).sum(1)==288) and np.all((A*B*g).sum(1)==0);assert len(np.unique(np.c_[A,B],axis=0))==8358
print('18D seed 8358 rational+sqrt3 weighted; verifying source',flush=True)
def bad(p,q):return ((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q))
for st in range(0,len(A),64):
 P=(A[st:st+64]*g)@A.T+3*(B[st:st+64]*g)@B.T-144;Q=(A[st:st+64]*g)@B.T+(B[st:st+64]*g)@A.T;v=bad(P,Q);ix=np.arange(st,min(st+64,len(A)));v[np.arange(len(ix)),ix]=False;assert not v.any(),('invalid',st)
 if st%1024==0:print('verified',st,flush=True)
np.savez_compressed('research/generic_reflection_campaign_20261105/cache/seed18.npz',A=A,B=B,g=g,norm=288);print('EXACT VERIFIED 18D 8358',flush=True)
