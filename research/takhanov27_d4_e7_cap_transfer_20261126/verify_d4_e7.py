import numpy as np, itertools, sys, collections, hashlib
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
D=Path(sys.argv[1]).resolve()
if (D/'common.py').is_file() and (D/'data/d27_certificate.npz').is_file():
    # Native source layout: k-nic/KissingNumbers@77f0860b/version1/certificates.
    cf=(D/'common.py').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(cf)).encode()+bytes([0])+cf).hexdigest()=='d1ae49a0f3ab6e9fb9906a98592ab9dbc6d9d7d5'
    sys.path.insert(0,str(D))
    from common import leech_shell
    fp=D/'data/d27_certificate.npz'
    assert hashlib.sha256(fp.read_bytes()).hexdigest()=='8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94'
    f=np.load(fp,allow_pickle=False);L=leech_shell().astype('<i2')
    keep=np.ones(len(L),dtype=bool);keep[f['deleted_leech_indices']]=False
    Q=np.vstack([np.pad(L[keep],((0,0),(0,3))),f['added_integer'].astype('<i2')])
else:
    for filename,digest in {'d27_certificate.npz':'8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94','L.npy':'68b3f12fa1d7024c6e3fafd91b3ea616526ccb1299b76c3b6e6e7d53d89ff302','Q.npy':'9674e0dc67eebf83be574c7da07292cd98a8a6c533b4d786591bb08390165d3b'}.items():
        assert hashlib.sha256((D/filename).read_bytes()).hexdigest()==digest
    f=np.load(D/'d27_certificate.npz');L=np.load(D/'L.npy');Q=np.load(D/'Q.npy')
assert np.array_equal(2*f['added_integer'][f['midpoint_rows'],:24],L[f['midpoint_endpoints']].sum(1))
assert len(L)==196560 and L.shape[1]==24
assert len(Q)==210108 and np.array_equal(Q[-26796:],f['added_integer']) and np.all(np.sum(Q.astype('i4')**2,1)==128)

# Pair up the 16 antipodal octagons by their Weyl-closed D4 carriers.
X=np.asarray(f['added_integer'],dtype=np.int16); E0=len(Q)-len(X);I=np.asarray(f['midpoint_rows'],int);ends=np.asarray(f['midpoint_endpoints'],int)
assert E0==183312 and len(I)==26496
def tup(r):return tuple(map(int,r))
edges=set(map(tuple,np.sort(ends,axis=1).tolist()));adj=collections.defaultdict(list)
for a,b in sorted(edges):adj[a].append(b);adj[b].append(a)
assert len(adj)==13248 and all(len(v)==2 for v in adj.values())
seen=set();cycles=[]
for a in sorted(adj):
    if a in seen:continue
    q=[];cur=a;prev=-1
    while cur not in seen:
        seen.add(cur);q.append(cur);nxt=next(v for v in adj[cur] if v!=prev);prev,cur=cur,nxt
    if len(q)==8 and np.array_equal(L[q[0]],-L[q[4]]):cycles.append(q)
assert len(cycles)==16
roots_in=set(map(tup,L))
def weyl(q):
    basis=L[q[:4]].astype(int);R=set(map(tup,basis));todo=list(R)
    for rr in todo:
        r=np.array(rr)
        for b in basis:
            dot=int(r@b);assert dot%64==0
            z=tup(r-(dot//64)*b)
            if z not in R:R.add(z);todo.append(z)
    assert len(R)==24 and R<=roots_in
    return tuple(sorted(R))
carriers=collections.Counter(weyl(q) for q in cycles);assert len(carriers)==8 and set(carriers.values())=={2}
H=set();frames=[]
for roots in sorted(carriers):
    R=np.array(roots,'i2');short=set()
    for a,b in itertools.combinations(R,2):
        if a@b==0:
            assert np.all((a+b)%2==0)
            h=(a+b)//2;assert int(h@h)==64;short.add(tup(h))
    assert len(short)==24
    H.update(short);F=[]
    for h in sorted(short):
        if all(np.dot(v,h)==0 for v in F):F.append(h)
    assert len(F)==4 and np.array_equal(np.array(F)@np.array(F).T,64*np.eye(4,dtype=int))
    frames.append(np.array(F,np.int16))
assert len(H)==160
# For any O(7) rotation of full E7 in each common D4-head/3-tail span.
# Integer projections use magnification by 8: the seven coordinates are (q.F,8*q_tail).
def color8(P):
    G=(abs(P@P.T)==2048) # actual projections have dot +/-32, magnification squared =64
    n=len(P);assert n==24 and not G.diagonal().any()
    order=sorted(range(n),key=lambda i:-int(G[i].sum()))
    C=[[] for _ in range(8)]
    def paint(k):
        if k==n:return True
        i=order[k]
        for c in C:
            if all(not G[i,j] for j in c):
                empty=not c;c.append(i)
                if paint(k+1):return True
                c.pop()
                if empty:break
        return False
    assert paint(0) and all(not G[np.ix_(c,c)].any() for c in C)
    return C
for ci,F in enumerate(frames):
    P=np.column_stack([Q[:,:24]@F.T,8*Q[:,24:]]).astype(np.int32);sq=np.einsum('ij,ij->i',P,P);assert sum(sq==8192)==68 and sum(sq==6144)==96
    lines=collections.Counter()
    for p in P[sq==6144]:
        sig=np.sign(p[np.flatnonzero(p)[0]]);lines[tup(sig*p)]+=1
    assert len(lines)==24 and set(lines.values())=={4}
    pp=np.array(sorted(lines),dtype=np.int32);colors=color8(pp)
    print('D4 carrier',ci,'colors on lexicographically sorted projected lines:',colors,flush=True)
# E7 Cartan normalization, squared norms of fundamental weights and highest root.
A=2*np.eye(7,dtype=int)
for i,j in [(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)]:A[i,j]=A[j,i]=-1
W=set(map(tup,np.eye(7,dtype=int)));todo=list(W)
for u in todo:
    q=np.array(u)
    for i in range(7):
        p=q.copy();p[i]-=A[i]@q;z=tup(p)
        if z not in W:W.add(z);todo.append(z)
m=np.array([2,3,4,3,2,1,2]);assert len(W)==126 and tup(m) in W and np.all(np.array(sorted(W))<=m)
from fractions import Fraction as Fr
def determinant(a):
    if len(a)==0:return 1
    return sum(((-1)**j)*x*determinant([r[:j]+r[j+1:] for r in a[1:]]) for j,x in enumerate(a[0]) if x)
A=A.tolist();det=determinant(A);D=[Fr(determinant([r[:j]+r[j+1:] for k,r in enumerate(A) if k!=j]),det) for j in range(7)]
assert D==[2,6,12,Fr(15,2),4,Fr(3,2),Fr(7,2)] and max(D[i]/int(m[i])**2 for i in range(7) if i!=5)<Fr(3,2)
# Strict-cap guards: every nonaxial point on one of the 160 short heads will pick one of twelve colors.
T=np.array(sorted(set(map(tup,X[I,24:]))),dtype='i2');assert len(T)==12 and np.all(np.sum(T**2,axis=1)==32)
parts=[I[np.all(X[I,24:]==t,axis=1)] for t in T]
tm={tup(t):k for k,t in enumerate(T)}
H=np.array(sorted(H),dtype='i2');D=H@X[I,:24].T
# Dots are exact int16 (input bounds <2^15). Convert actual X row numbers to position in D.
old_to_mid={int(i):j for j,i in enumerate(I)}
idx=[np.array([old_to_mid[int(j)] for j in l],int) for l in parts]
G=[]
for h in range(len(H)):
    # Find one pair H0 of coordinate axes such that all four signed channels on that pair have maximum >=64.
    choices=[d for d in range(3) if all(D[h,idx[k]].max()>=64 for k,t in enumerate(T) if t[d]==0)]
    assert choices; d=choices[0];pair=set(range(3))-{d}
    for k,t in enumerate(T):
        guards=set(map(int,idx[k][D[h,idx[k]]>=32]))
        if t[d]!=0:
            i=next(j for j in pair if t[j]!=0);j=next(j for j in pair if t[j]==0)
            aux=np.zeros(3,dtype='i2');aux[i]=t[i];aux[j]=t[i] if i<j else -t[i];ell=tm[tup(aux)]
            guards.update(map(int,idx[ell][D[h,idx[ell]]>=64]))
        G.append(guards)
assert len(G)==1920
delidx=list(map(int,f['deleted_leech_indices']));di={r:k for k,r in enumerate(delidx)}
refunds=[set() for r in delidx]
for i,(a,b) in enumerate(ends):
    refunds[di[int(a)]].add(i);refunds[di[int(b)]].add(i)
assert len(refunds)==13248 and all(len(r)==4 for r in refunds)
G+=refunds;ind=[0];flat=[]
for g in G:flat.extend(sorted(g));ind.append(len(flat))
mat=csr_matrix((np.ones(len(flat),dtype=np.int8),flat,ind),shape=(len(G),len(I)))
M=maximum_bipartite_matching(mat,perm_type='column')
assert len(M)==15168 and len(set(map(int,M)))==15168 and min(M)>=0 and all(int(j) in g for j,g in zip(M,G))
# Disjoint source slices for optional arbitrary-headed axial and pure-tail trades.
for j in range(3):
    for s in [-8,8]:
        t=np.zeros(3,dtype='i2');t[j]=s
        Z=X[np.all(X[:,24:]==t,axis=1),:24];R=Z@Z.T
        assert len(Z)==48 and np.all(np.diag(R)==64) and np.all(np.sum(R==-64,axis=1)==1) and set(np.unique(R))=={-64,0,64}
Z=X[np.all(X[:,:24]==0,axis=1)];assert len(Z)==12 and np.all(np.sum(Z[:,24:]**2,axis=1)==128)
print('ALL EIGHT: >=132 source blockers for 126 rotated E7 roots. Continuous short-head code: 15168 distinct native midpoint blockers assigned; no gain.')

# If the new body lies farther out on any short-head ray (arbitrary radial factors), 
# up to four new sites on that ray can be charged to four distinct retained equatorial roots.
keep=np.ones(len(L),dtype=bool);keep[f['deleted_leech_indices']]=False
RG=[];allct=[]
for h in H:
    v=L@h;assert np.max(v)==64 and np.sum(v==64)==46
    r=list(map(int,np.flatnonzero((v==64)&keep)));allct.append(len(r));RG.extend([r]*4)
assert collections.Counter(allct)=={38:152,22:8}
ind=[0];jj=[]
for r in RG:jj.extend(r);ind.append(len(jj))
Aeq=csr_matrix((np.ones(len(jj),dtype=np.int8),jj,ind),shape=(len(RG),len(L)))
MM=maximum_bipartite_matching(Aeq,perm_type='column')
assert len(MM)==640 and min(MM)>=0 and len(set(map(int,MM)))==640 and all(int(j) in g for j,g in zip(MM,RG))
print('Arbitrary outward displacement: 640 different retained equatorial blockers for 4 slots on each of 160 rays.')
