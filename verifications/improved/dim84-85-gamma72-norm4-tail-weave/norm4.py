"""Finite geometry check for the added norm-four tails."""
from pathlib import Path
import json
from itertools import combinations,product
import numpy as np
def norm4(n,K,P):
    M,H,R,tail,pos=(K[f]for f in ('M','H','R','tail','pos'))
    H=list(map(int,H))
    tail(n)
    PART=((5,0),(9,6),(3,2),(11,7),(4,1),(10,8)) if n==84 else ((10,6),(3,11),(9,4),(0,7),(5,2),(1,8))
    assert sorted(sum(PART,()))==list(range(12))
    GROUP=[PART[2*x]+PART[2*x+1]for x in range(3)]
    Y={x:tuple(sorted(y for j,Q in enumerate(PART)if (j^i)&1 and j//2!=i//2 for y in Q)) for i,Q in enumerate(PART) for x in Q}
    edges={tuple(sorted((x,y)))for x in Y for y in Y[x]}
    assert len(edges)==24 and all(len(Y[x])==4 and all(x in Y[y]for y in Y[x]) for x in Y)
    assert all(not(set(Y[x])&set(Y[y]))for x,y in edges)
    W=[]; phase=[]
    def put(cluster):
        phase_id=max(phase,default=0)+1
        for v in cluster:
            z=np.zeros(n-72,dtype='i8')
            for i,e in v: z[i]=e
            W.append(z);phase.append(phase_id)
    # 168 norm-four D12 roots: delete the four on each of 24 paired edges.
    for i in range(11):
        mp=[(11,i)]+[((i+j)%11,(i-j)%11)for j in range(1,6)]
        assert sorted(v for p in mp for v in p)==list(range(12))
        put([((a,2*u),(b,2*v))for a,b in mp if tuple(sorted((a,b)))not in edges for u,v in product((-1,1),repeat=2)])
    assert len(W)==168
    # 384 octad tails: only even sign words, split into affine RM(1,3) cosets.
    RM={sum((((i&a).bit_count()^b)&1)<<i for i in range(8))for a in range(8)for b in range(2)}
    assert len(RM)==16 and {w.bit_count() for w in RM}=={0,4,8}
    for A,B in combinations(GROUP,2):
        s=A+B; U={i for i in range(256)if i.bit_count()%2==0}
        while U:
            w=min(U); C={w^v for v in RM};assert C<=U;U-=C
            put([tuple((i,1-2*((w>>j)&1))for j,i in enumerate(s))for w in sorted(C)])
    assert len(W)==552
    # 192 2+1+1+1+1 tails. For each major coordinate x, group both major
    # signs together with two complementary even signs on the neighbor quartet.
    for x in sorted(Y):
        for mask in range(8):
            if mask.bit_count()%2:continue
            put([((x,e),)+tuple((y,1-2*((w>>j)&1))for j,y in enumerate(Y[x])) for e in (-2,2) for w in (mask,mask^15)])
    assert len(W)==744
    if n==85:
        # On each disjoint coordinate quartet, all 16 signings of tail
        # coordinates consisting of the 13th doubled and four signed units.
        for group in GROUP:
            for mask in range(8):
                if mask.bit_count()%2:continue
                put([((12,e),)+tuple((y,1-2*((w>>j)&1))for j,y in enumerate(group)) for e in (-2,2) for w in (mask,mask^15)])
    X=np.array(W,dtype='i8'); T=np.array(phase)
    count,phases=(744,83)if n==84 else (792,95)
    assert X.shape==(count,n-72) and T.min()==1 and T.max()==phases and len(set(T))==phases
    assert len(np.unique(X,axis=0))==len(X)
    if n==84:
        J=np.array(json.loads((P/'tail84.json').read_text()),dtype='i8')
        assert J.shape==(420,12)
        assert np.max(np.abs(X@J.T))<=4
    else:
        def parse(t):
            return (0,int(t.split('*')[0]))if '*' in t else (int(t),0)
        J=np.array([[parse(t)for t in r.split(',')]for r in (P/'ze99.txt').read_text().splitlines()],dtype='i8')
        assert J.shape==(1154,13,2) and np.max(abs(J))<=4
        A=X@J[:,:,0].T; B=X@J[:,:,1].T
        assert not (pos(A-8,B,3)|pos(-A-8,-B,3)).any()
    # Both tests give physical old mixed tail dot <=2; pure tail dot <=4.
    G=X@X.T
    assert np.array_equal(np.diag(G),np.full(len(G),8))
    G-=8*np.eye(len(G),dtype='i8')
    assert G[T[:,None]==T[None,:]].max()<=0
    assert G[T[:,None]!=T[None,:]].max()<=4
    # Actual horizontal shell for every used B-only nonzero phase: 24 antipodal
    # raw Leech norm8 vectors in each of three blocks, physical norm4.
    D=dict(zip(map(tuple,R["A8"]),R["V8"]))
    for h in T:
        V=D[0,H[h]]
        assert V.shape==(24,24) and np.array_equal(V@R['F']@V.T,8*np.eye(24,dtype='i8'))
    assert all(H[i]!=0 for i in T)
    # B-only difference clears distinct nonzero phase pairs by horizontal d²>=4.
    # Difference from every old mixed color is a nonzero first color at H,
    # clearing by d²>=6. Same new phase but different rows clears vertically.
    return 144*len(X)
