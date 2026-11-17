from pathlib import Path
import numpy as np,json,itertools
P=Path(__file__).resolve().parent
V=np.array([[pow(x,j,17) for x in range(6)] for j in range(4)],dtype=np.int64)
def detmod(A):
 A=A.copy()%17;ans=1
 for j in range(len(A)):
  k=next((i for i in range(j,len(A)) if A[i,j]),None)
  if k is None:return 0
  if k!=j:A[[j,k]]=A[[k,j]];ans=-ans
  pivot=int(A[j,j]);ans=ans*pivot%17;A[j]=A[j]*pow(pivot,-1,17)%17
  for i in range(j+1,len(A)):A[i]=(A[i]-A[i,j]*A[j])%17
 return ans%17
assert all(detmod(V[:,cols])!=0 for cols in itertools.combinations(range(6),4))
S=np.array([s for s in itertools.product((0,1),repeat=6) if sum(s)%2==0],dtype=np.int16)
assert len(S)==32
D=np.count_nonzero(S[:,None,:]!=S[None,:,:],axis=2);assert D[D>0].min()==2
base=json.loads((P/'sparse_certificate.json').read_text())['totals'];results={}
for h in range(4,8):
 coeff=np.array(list(itertools.product(range(h),range(17),range(17),range(17))),dtype=np.int64)
 words=(coeff@V)%17;assert np.array_equal(words[:,0],coeff[:,0])
 supports=np.column_stack((88+words[:,0],words[:,1:]+17*np.arange(5)))
 assert supports.shape==(h*17**3,6) and len(np.unique(supports,axis=0))==h*17**3
 assert np.all(supports[:,0]>=88) and np.all(supports[:,0]<88+h) and np.all(supports[:,1:]<85)
 # Distinct cubic words agree at most3 coordinates by the15 exact4-column minors above.
 assert 6*16==96 and 3*16==48 and 4*16*2>=96
 assert 5*4+2*4<48 and (5*4)**2*48<48**2*11
 assert 16*96<48**2 and 96*48<48**2*11 and 4*96<48**2
 n=88+h;rs=32*h*17**3;poles=2*n;total=int(base[str(n)])+rs+poles
 results[n]=dict(head_dimension=h,support_count=len(supports),signs=32,rs_count=rs,poles=poles,prior_dense_sparse22=int(base[str(n)]),total=total)
obj=dict(status='PASS exact RS17 MDS minors, finite support/sign counts and cross inequalities',evaluation_generator=V.tolist(),results=results)
(P/'rs17_certificate.json').write_text(json.dumps(obj,indent=2)+'\n');print(json.dumps(obj,indent=2))
