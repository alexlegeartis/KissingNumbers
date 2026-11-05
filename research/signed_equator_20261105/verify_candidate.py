import numpy as np,sys
z=np.load(sys.argv[1]);den=int(z['den']) if 'den' in z else 1;A=z['rational'].astype(np.int64);B=z['radical'].astype(np.int64);assert A.shape==B.shape and A.shape[1]==13 and len(A)>1154
p=(A*A+3*B*B).sum(1);q=2*(A*B).sum(1);assert np.all(p==16*den*den) and np.all(q==0)
n=len(A);bad=0
for i in range(0,n,128):
 p=A[i:i+128]@A.T+3*B[i:i+128]@B.T-8*den*den;q=A[i:i+128]@B.T+B[i:i+128]@A.T;bb=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q));rr=np.arange(i,i+len(p))[:,None];cc=np.arange(n)[None,:];bad+=np.count_nonzero(bb&(rr<cc))
print('points',n,'bad',bad,flush=True);assert bad==0;print('CERTIFIED 13D >=',n)
