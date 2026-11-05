"""Check full 4C orbit expansion and target distance cuts by direct arithmetic."""
import numpy as np
G=np.zeros((12,24),np.int32);g=np.array([7,1,4,4,2,7,3,3,4,7,2,1])
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
gray=np.array([sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)],np.uint8);wt=np.array([int(i).bit_count() for i in range(16)])
t=np.arange(4096,dtype=np.uint16);bits=np.stack([(t>>i)&1 for i in range(12)],axis=1).astype(np.int32);golay=(bits@G)%2;assert len(np.unique(golay,axis=0))==4096 and np.min(golay[1:].sum(1))==8
b=np.array(list('305322310155'),np.int32);source=gray[(b@G)%8];expanded=source[None,:]^(15*golay).astype(np.uint8);direct=gray[((b[None,:]+4*bits)@G)%8];assert np.array_equal(expanded,direct) and np.all(wt[expanded[:,:4]]==2) and len(np.unique(expanded,axis=0))==4096
target=np.array([int(x,16) for x in '6596e71dd17ebbb71212bee8'],np.uint8);assert np.all(wt[target[:4]]==2) and np.all(wt[target[4:]]%2==1);dist=wt[expanded^target].sum(1);assert np.all(dist>=20) and np.all(dist%2==0)
print('PASS 4096 exact source translates, distance cuts and candidate orbit min distance 32; no compatibility proof')
