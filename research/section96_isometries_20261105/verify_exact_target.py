"""Independent vectorized exhaustive check of one mixed-parity target."""
import numpy as np,sys
word=sys.argv[1];target=np.array([int(x,16) for x in word],np.uint8);assert len(target)==24
gray=np.array([sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)],np.uint8);wt=np.array([int(x).bit_count() for x in range(16)]);assert np.all(wt[target[:4]]==2) and np.all(wt[target[4:]]%2==1)
g=np.array([7,1,4,4,2,7,3,3,4,7,2,1]);G=np.zeros((12,24),np.int32)
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
choices=[]
for j in range(12):
 d=wt[gray^target[j]]
 if j<4:
  choices.append((np.flatnonzero((wt[gray]==2)&(d==0)),np.flatnonzero((wt[gray]==2)&(d==2))))
 else:choices.append((np.flatnonzero(d==1),np.flatnonzero(d==3)))
assert all(len(choices[j][0])==1 and len(choices[j][1])==4 for j in range(4));assert all(len(choices[j][0])==len(choices[j][1])==4 for j in range(4,12))
t=np.arange(4**8,dtype=np.uint32);leaves=0;bad=0;first=None
for pos in range(-1,12):
 for replacement in range(4 if 0<=pos<4 else 1):
  z=np.zeros((len(t),12),np.int32)
  for j in range(4):z[:,j]=choices[j][1][replacement] if pos==j else choices[j][0][0]
  for j in range(4,12):z[:,j]=choices[j][1 if pos==j else 0][(t>>(2*(j-4)))&3]
  a=np.zeros_like(z)
  for j in range(12):a[:,j]=(7*(z[:,j]-a[:,:j]@G[:j,j]))%8
  tail=gray[(a@G[:,12:])%8];dist=wt[tail^target[12:]].sum(1);firstcost=wt[gray[z]^target[:12]].sum(1);total=firstcost+dist;ix=np.flatnonzero(total<24);bad+=len(ix);leaves+=len(t)
  if len(ix) and first is None:first=(a[ix[0]].tolist(),int(total[ix[0]]))
assert leaves==1638400
print('leaves',leaves,'conflicts',bad,'first',first,'EXACT_HOLE' if not bad else 'COVERED')
