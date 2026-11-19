"""Independent signed 38D triangle construction and changed-pair geometry."""
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
import numpy as np, json
D=Path(__file__).resolve().parent
P=D.parent/'dim38-joint-latitudes'
S=P/'source/source'
# Construct the four triangles as literal signed D5 root arrays. Each row
# gives three oriented roots, scaled by four only at embedding time.
T=[
 [[1,-1,0,0,0],[0,1,-1,0,0],[-1,0,1,0,0]],
 [[1,1,0,0,0],[0,-1,0,-1,0],[-1,0,0,1,0]],
 [[1,0,1,0,0],[0,0,-1,0,-1],[-1,0,0,0,1]],
 [[0,1,0,-1,0],[0,0,0,1,1],[0,-1,0,0,-1]],
]
T=np.array(T,dtype=np.int64)
assert all(np.array_equal(x@x.T,3*np.eye(3,dtype=np.int64)-1) for x in T)
X=[];Y=[]
for bits in product([0,1],repeat=7):
    syndrome=sum((sum(i for j,i in enumerate(bits,1) if j>>k&1)&1)<<k for k in range(3))
    p=sum(bits)&1
    b=14+5*p
    for t in T[syndrome>>1]:
        body=np.zeros(24,dtype=np.int64);body[b:b+5]=4*t*(-1 if syndrome&1 else 1)
        word=1-2*np.array(bits,dtype=np.int64)
        for tail in (np.r_[word,np.zeros(7,dtype=np.int64)],np.r_[np.zeros(7,dtype=np.int64),word]):
            X.append(body);Y.append(tail)
B=np.array(X);H=np.array(Y)
assert len(B)==len(np.unique(np.c_[B,H],axis=0))==768
# Exact paired dot numerators: unit dot equals B.B'/96 + (2/21) H.H'.
g=B@B.T;z=H@H.T
assert np.array_equal(np.diag(g),np.full(768,32)) and np.array_equal(np.diag(z),np.full(768,7))
assert (np.triu(21*g+192*z,1)<=1008).all()
# Independent enumeration of Golay/Leech line indexing.
G=np.array([[int(c) for c in s] for s in (S/'golay_basis.txt').read_text().splitlines()],dtype=np.int64)
C=[np.zeros(24,dtype=np.int64)]
for row in G:C.extend([x^row for x in C])
A=[]
for i,j in combinations(range(24),2):
 for sg in [-4,4]:
    a=[0]*24;a[i]=4;a[j]=sg;A.append(a)
for word in C:
 I=np.flatnonzero(word)
 if len(I)==8:
    for signs in product([-2,2],repeat=7):
        if sum(x<0 for x in signs)%2:continue
        a=[0]*24
        for j,v in zip(I,(2,*signs)):a[j]=v
        A.append(a)
 sg=1-2*word
 for i in range(24):
    a=sg.copy();a[i]*=-3
    if a[0]>0:A.append(a.tolist())
A=np.array(sorted(A),dtype=np.int64)
assert A.shape==(98280,24)
assert np.all(np.sum(A*A,axis=1)==32)
L={tuple(a):i for i,a in enumerate(A)}
R=set()
for a in B:
 i=int(np.flatnonzero(a)[0]);R.add(L[tuple(a if a[i]>0 else -a)])
assert len(R)==24
old=json.loads((P/'refill/certificate.json').read_text())
prev=json.loads((P/'hamming/certificate.json').read_text())
gd=np.array(old['body_line_indices'],dtype=np.int64)
ng=set(gd.tolist());assert len(ng)==14 and ng.isdisjoint(R)
C0=A[[i for i in range(len(A)) if i not in R|ng]]
equator=np.loadtxt(S/'equator_lines_6855.txt',dtype=np.int64)
W=np.loadtxt(S/'dim14_directions.txt',delimiter=',',dtype=np.int64)
# One independently checked rational upper envelope on the positive
# coefficients of all cap interactions, with kg=48 sqrt(2), kh=6 sqrt(7).
assert 67**2<48**2*2 and F(79,5)**2<36*7
ub=int(np.abs(np.unique(B,axis=0)@C0.T).max())
uh=int(np.abs(H@W.T).max())
assert F(ub,67)+F(5*uh,79)<F(1,2)
# Compare squares where one orthogonal block suffices.
grid=np.array(old['doubled_head_lines'],dtype=np.int64)*2
axis=W[old['axis_keep']]@np.array(old['rotation_numerator'],dtype=np.int64).T
eq=int(np.abs(B@equator.T).max())
assert np.all((B@A[gd].T)==0)
gmax=int(np.abs(H@grid.T).max())
ax=int(np.abs(H@axis.T).max())
assert F(eq*eq,96*48)<=F(1,4)
assert F(2*gmax*gmax,21*128)<=F(1,4)
assert F(2*ax*ax,21*8)<=F(1,4)
assert F(2,21)<=F(1,4)
rest=sorted(set(prev['removed_cap_body_lines'])-(R|ng))
assert len(rest)==5 and np.all((A[rest]@A[gd].T)==0)
count=6*len(C0)+2*len(equator)+len(axis)+28+4*len(grid)+len(B)
assert count==611798
print(json.dumps({'status':'PASS','dimension':38,'total':count,'new_points':len(B),'new_pairs':len(B)*(len(B)-1)//2,
 'maxima':{'cap_body':ub,'cap_head':uh,'equator':eq,'grid_tail':gmax,'axis_tail':ax}}))
