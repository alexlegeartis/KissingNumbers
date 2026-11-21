"""Independent replay of layered syndrome rescue. No original verifier imported."""
from contextlib import redirect_stdout
from io import StringIO
from itertools import combinations
from collections import Counter, defaultdict
with redirect_stdout(StringIO()):
    from audit import G,D,sign,receipt,data,X,Y,dot
wt=int.bit_count
blockset=[[i for i in range(12) if i//4!=t] for t in range(3)]
def proj(w,J):
    return sum(((w>>(2*i))&3)<<(2*j) for j,i in enumerate(J))
perms=[]
for t,J in enumerate(blockset):
    C={proj(w,J) for w in D}
    assert len(C)==4096
    # Highest pivots pick the numerically smallest representative of each coset.
    # Find the 16 representatives and a direct 65536-entry quotient table
    # without using the original elimination.
    q=[-1]*65536
    reps=[]
    for v in range(65536):
        if q[v]!=-1: continue
        reps.append(v)
        for c in C:
            assert q[v^c]==-1
            q[v^c]=v
    assert len(reps)==16
    holes=[i for i in range(16) if (1<<i) in reps]
    encode=lambda a:sum(((a>>j)&1)<<i for j,i in enumerate(holes))
    ids=[0,5,9,12] if t==1 else [0,5,10,15]
    v=[encode(a) for a in ids]
    assert all(x in reps for x in v)
    for x,y in combinations(v,2):
        assert min(sum(bool((w>>(2*j))&3) for j in range(8)) for c in C for w in [x^y^c])>=2
    perms.append((J,q,v))
def syn(w):
    return tuple(q[proj(w,J)] for J,q,_ in perms)
L={syn(0):0}
for i in range(24):
    L.update({syn(x^(1<<i)): x^(1<<i) for x in list(L.values())})
assert len(L)==4096 and all(syn(d)==(0,0,0) for d in D)
O=defaultdict(list)
for x in G:
    if wt(x)==8:
        for s in combinations([j for j in range(24) if (x>>j)&1],4):
            O[sum(1<<j for j in s)].append(x)
assert len(O)==10626 and set(map(len,O.values()))=={5}
GS=set(G)
assert all(abs(dot(a,b))<=4 for T in X for a in T for S in Y for b in S)
prior={92:30151281256,93:44038752920,94:74103655116,95:110610877508,96:209198254912}
for line in data('shifted_patterns.json'):
    mask=line['mask']
    ids=[int(x,36) for x in line['patterns'].split(',')]
    assert len(set(ids))==len(ids) and all(0<=x<2**24 for x in ids)
    Q=[]
    for idx in ids:
        a,b=G[idx%4096],G[idx//4096]
        assert wt(a|b)==12 and not (mask&(a|b))
        T=(a&~b,b&~a,a&b)
        sizes=sorted(map(wt,T))
        assert sizes in ([4,4,4],[0,0,12])
        if sizes==[4,4,4]:
            q=[]
            for i in range(3):
                x,y=T[(i+1)%3],T[(i+2)%3]
                F=sorted(c^x for c in O[x] if c&(x|y)==x)
                assert len(set(F))==4 and T[i] in F and all(f|y in GS for f in F)
                assert all(not x&y for x,y in combinations(F,2))
                q.append(perms[i][2][F.index(T[i])])
            S=L[tuple(q)]
        else: S=0
        Q.append((a,b,T,S))
    hist=Counter()
    rescued=[0,0,0]
    for (a,b,T,A),(c,d,V,B) in combinations(Q,2):
        dist=wt((a^c)|(b^d)); hist[dist]+=1
        if dist>=12: continue
        assert dist==8 and all(wt(x)==4 for x in (*T,*V))
        J=[i for i in range(3) if T[i]!=V[i]]
        assert len(J)==1
        u=J[0]
        assert not T[u]&V[u] and wt((a|b)&(c|d))==8
        J,q,v=perms[u]
        p=q[proj(A^B,J)]
        assert p in {x^y for x,y in combinations(v,2)}
        rescued[u]+=1
    assert len(D)==4096 and len(sign)==256
    print(mask,len(Q),sum(rescued),rescued,dict(sorted(hist.items())))
    layer=len(Q)*len(D)*len(sign)
    for n in ({15:[92],3:[93,94,95],0:[96]}[mask]):
        total=receipt[n]['incumbent']+layer
        print(n, layer, total, total-prior[n])
