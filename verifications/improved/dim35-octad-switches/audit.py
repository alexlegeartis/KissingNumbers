"""Independent exact 35D audit by explicitly expanding optional octad parities."""
import json,itertools as I, numpy as np
from pathlib import Path
if not __debug__:
    raise RuntimeError("Run without -O")
P=Path(__file__).resolve().parent
Q=P.parent/'dim35-shared-root-heads'
J=json.loads((P/'witness.json').read_text());C=json.loads((Q/'code_data.json').read_text())
Z=np.load(Q/'shared_positive.npz',allow_pickle=False)
old=Z['new_stars'].astype(np.int64);H=Z['heads'].astype(np.int64);D=Z['dense_positions'].astype(np.int64)
assert sorted(D.tolist())==[i for i in range(35)if i not in (0,3,18)]
S=[tuple(i for i in range(35)if int(x)>>i&1) for x,_ in C['full_parity_supports']]
used=set()
for i,u,v in J['octads']:
 assert tuple(u)==S[i] and i not in used;used.add(i);S[i]=tuple(v)
assert len(S)==len(set(map(frozenset,S)))==2202 and all(len(set(s))==8 and 0<=min(s)<=max(s)<35 for s in S)
assert max(len(set(s)&set(t)) for s,t in I.combinations(S,2))<=4
r=J['remove'];a,b=r['start'],r['stop']
cx,*ws=r['template'];assert set(np.unique(np.argwhere(old[a:b])[:,1]))=={cx,*ws}
assert np.all(old[a:b,cx]==4)
assert all(set(np.flatnonzero(x))=={cx,*ws} for x in old[a:b])
pairs=[]
def star(u,c):
 a,*L=u
 for row in I.product((-2,2),repeat=4):
  if sum(q<0 for q in row)%2:continue
  v=np.zeros(35,dtype=np.int64);v[a]=c;v[L]=row;pairs.append(v)
for t in J['both']:
 for q in (-4,4):star(t,q)
star(J['positive'],4)
T=np.vstack([old[:a],old[b:],pairs]); assert len(T)==246
assert np.all(np.sort(abs(T),axis=1)==[0]*30+[2]*4+[4])
assert np.all(np.sort(abs(H),axis=1)==[1]*33+[3]*2)
G=np.r_[T,H];N=np.sum(G*G,axis=1)
assert len(set(map(tuple,G)))==len(G)
def check(A,B,n,m,diag=False):
 q=A@B.T
 if diag:np.fill_diagonal(q,0)
 q=np.maximum(q,0)
 assert np.all(4*q*q<=np.asarray(n)*np.asarray(m))
check(G,G,N[:,None],N[None,:],True)
R=[]
for x,y in I.combinations(range(35),2):
 for u,v in I.product((-4,4),repeat=2):
  z=np.zeros(35,dtype=np.int64);z[[x,y]]=[u,v];R.append(z)
R=np.array(R);assert len(R)==2380
M=G@R.T;bad=np.any(4*np.maximum(M,0)**2>32*N[:,None],axis=0);R=R[~bad];assert len(R)==2172
# Enumerate every signing on every octad. Choose a compatible parity.
Pats=np.array(list(I.product((-2,2),repeat=8)),dtype=np.int64)
P0=(Pats<0).sum(1)%2
prior=dict(json.loads((Q/'full_parities.json').read_text()))
X=G
for i,s in enumerate(S):
 w=X[:,s] @ Pats.T
 badp=np.any(4*np.maximum(w,0)**2>32*N[:,None],axis=0)
 good=[p for p in (0,1) if not np.any(badp[P0==p])]
 assert good, s
 m=sum(1<<j for j in s)
 p=0 if i in used else prior[m]
 assert p in good
# Explicit code expansion with exact Hamming-distance verification and
# all star/head products. Ordinary octads and roots pass by magnitude.
W=[0]
assert all(0<int(t)<(1<<32) for t in C["dense_generator"]) and 0<=int(C["dense_coset_shift"])<1<<32
for t in C['dense_generator']:W.extend([int(t)^x for x in W])
assert len(W)==len(set(W))==131072 and min(w.bit_count()for w in W if w)>=8
A=np.asarray(W,dtype=np.int64)^int(C['dense_coset_shift'])
B=1-2*((A[:,None]>>np.arange(32))&1)
for block in np.array_split(B,32):
 Y=X[:,D]@block.T
 assert np.all(4*np.maximum(Y,0)**2<=32*N[:,None])
total=len(B)+128*len(S)+len(T)+len(H)+len(R)
assert total==415358
print(json.dumps(dict(status='PASS',total=total,dense=len(B),octads=len(S),octad_points=128*len(S),stars=len(T),heads=len(H),roots=len(R),removed_roots=int(bad.sum())),sort_keys=True))
