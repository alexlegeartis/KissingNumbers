from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
G=np.zeros((8,14),dtype=np.int64)
g=[1,1,0,1,0,1,1]
for i in range(8):G[i,i:i+7]=g
# Independent finite proof: g divides x^14+1 over F3.
f=[1]+[0]*13+[1]
for i in range(14,5,-1):
 a=f[i]%3
 for j,b in enumerate(g):f[i-6+j]=(f[i-6+j]-a*b)%3
assert not any(f)
# The leading triangular 8x8 block proves row independence.
assert np.all(np.diag(G[:,:8])==1) and not np.tril(G[:,:8],-1).any()
D=np.array([[(i//3**j)%3 for j in range(8)] for i in range(3**8)],dtype=np.int64)
C=(D@G)%3
assert len(np.unique(C,axis=0))==6561
hist=np.bincount(np.count_nonzero(C,axis=1),minlength=15)
assert hist[0]==1 and not hist[1:5].any() and hist[5]>0
# Last coordinate sees only final generatorrow, so shortening simply drops that row and coordinate.
H=G[:7,:13];E=np.array([[(i//3**j)%3 for j in range(7)] for i in range(3**7)],dtype=np.int64);B=(E@H)%3
assert np.array_equal(B,C[C[:,-1]==0,:-1])
h13=np.bincount(np.count_nonzero(B,axis=1),minlength=14)
assert len(np.unique(B,axis=0))==2187 and h13[0]==1 and not h13[1:5].any() and h13[5]>0
for n,a in [('generator14',G),('generator13',H),('phase14',C),('phase13',B)]:np.save(P/(n+'.npy'),a.astype(np.uint8))
report=dict(status='PASS',sources=['https://www.codetables.de/BKLC/BKLC.php?q=3&n=14&k=8','https://www.codetables.de/BKLC/BKLC.php?q=3&n=13&k=7'],polynomial_coefficients_ascending=g,weight_histogram14=hist.tolist(),weight_histogram13=h13.tolist(),phase14=6561,phase13=2187,addition91=2*31680*6561,addition92=4*31680*6561+4*6930*2187,total91=21943640934+2*31680*6561+2,total92=21943640934+4*31680*6561+4*6930*2187+6)
(P/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
