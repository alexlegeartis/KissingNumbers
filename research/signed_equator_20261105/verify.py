import numpy as np,collections
from pathlib import Path
z=np.load(Path(__file__).with_name('source.npz'));X=np.c_[z['rational'],z['radical']].astype(np.int64);A=X[:,:13];B=X[:,13:];assert X.shape==(1154,26)
# Verify the entire source in Q(sqrt(3)), including norms and all unordered pairs.
assert np.all((A*A+3*B*B).sum(1)==16) and np.all((A*B).sum(1)==0)
for st in range(0,len(X),128):
 p=A[st:st+128]@A.T+3*B[st:st+128]@B.T-8;q=A[st:st+128]@B.T+B[st:st+128]@A.T;bad=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q));ii=np.arange(st,st+len(p))[:,None];jj=np.arange(len(X))[None,:];assert not np.any(bad&(ii<jj))
z=(B==0).all(1);E=z&(A[:,12]==0)&(np.count_nonzero(A[:,:12],axis=1)==4)&np.isin(A[:,:12],[-2,0,2]).all(1);H=z&(np.abs(A[:,:12])==1).all(1)&(np.abs(A[:,12])==2);tail=z&(np.count_nonzero(A,axis=1)==1)&(np.abs(A[:,12])==4);assert (E.sum(),(H&(A[:,12]==2)).sum(),(H&(A[:,12]==-2)).sum(),tail.sum(),(~(E|H|tail)).sum())==(816,144,144,2,48)
# Every global sign word sees exactly 51 signed equator vectors.
words=np.arange(4096,dtype=np.uint16);bits=np.arange(12,dtype=np.uint16);signs=1-2*((words[:,None]>>bits)&1).astype(np.int16);U=A[E,:12];matches=((U[None,:,:]==0)|(U[None,:,:]*signs[:,None,:]>0)).all(2).sum(1);assert np.all(matches==51)
O=X[~(E|H|tail)];axes=np.zeros((24,13),np.int64)
for i in range(12):axes[2*i,i]=4;axes[2*i+1,i]=-4
p=axes@O[:,:13].T-8;q=axes@O[:,13:].T;bad=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q));assert np.all(bad.sum(1)==2) and np.all(bad.sum(0)==1)
extras=np.load(Path(__file__).with_name('extras.npy')).astype(np.int64);assert extras.shape==(69,26);assert np.all((extras[:,:13]**2+3*extras[:,13:]**2).sum(1)==16) and np.all((extras[:,:13]*extras[:,13:]).sum(1)==0)
expected=np.zeros((26,26),np.int64)
for i in range(13):expected[2*i,i]=4;expected[2*i+1,i]=-4
assert np.array_equal(extras[:26],expected)
more=np.load(Path(__file__).with_name('extras_den5.npy')).astype(np.int64);assert more.shape==(50,26);assert np.all((more[:,:13]**2+3*more[:,13:]**2).sum(1)==400) and np.all((more[:,:13]*more[:,13:]).sum(1)==0)
print('Verified: 1154-point source; 816 equator, 144+144 heads, 2 tail axes, 48 outsiders; 51-per-sign and 24 disjoint axis blocker pairs.')
