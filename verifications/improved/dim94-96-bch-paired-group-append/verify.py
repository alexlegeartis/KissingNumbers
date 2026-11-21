"""Independent exact replay for paired-group BCH support append to 94–96D."""
import itertools as it
import functools as fn
import operator as op
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def product(x,y):
    z=0
    for k in range(5):
        if (y>>k)&1:z^=x
        x*=2
        if x>=32:x^=37
    return z
assert all(len({product(a,b)for b in range(32)})==32 for a in range(1,32))
def make_H(x): return 1+2*x+64*product(x,product(x,x))
H=list(map(make_H,range(24)))
G=[H[k:k+4]for k in range(0,24,4)]
def subspace(L):
    S={0}
    for x in L:S|={s^x for s in S}
    return S
assert len(subspace(H))==2048
X=lambda L:fn.reduce(op.xor,L,0)
for k in range(1,6):assert all(X(L)!=0 for L in it.combinations(H,k))
F=[(0,0),(4,1025),(913,1987),(1759,327),(417,1632),(772,1281)]
def phase(i,a):
    p,q=F[i]
    return (p if a&1 else 0)^(q if a&2 else 0)
for i,j in it.combinations(range(6),2):
    gone=subspace(G[i]+G[j])
    survivors={0}|{b for k in range(6)if k not in (i,j)for b in G[k]}
    for a in (1,2,3):
        r=phase(i,a)^phase(j,a)
        assert all(r^q not in gone for q in survivors)
# Dense and older shifted rows have abs-dot<=4 with each flat row.
R=list(it.product((-1,1),repeat=4))
A=list(it.product((-1,1),repeat=4))
for i,j in it.combinations(range(4),2):
    A.extend(tuple(a if k==i else b if k==j else 0 for k in range(4)) for a,b in it.product((-2,2),repeat=2))
A.extend(tuple(a if j==i else 0 for j in range(4))for i in range(4)for a in (-2,2))
assert all(abs(sum(a*b for a,b in zip(x,y)))<=4 for x,y in it.product(R,A))
assert 4*(6*4)**2==24*96 and 4*48<24**2
baseline={94:74136160972,95:117085834308,96:209283189568}
records=[]
for dim in (94,95,96):
    sizes=[4 if dim==96 else 2]+[4]*5
    offsets=[sum(sizes[:k])for k in range(6)]
    # Virtual slot i chooses one block of physical group (i//2)*2+(a//2),
    # in its first or second half as prescribed by i mod 2.
    def row(i,a):
        group=2*(i//2)+a//2
        half=sizes[group]//2
        if a%2>=half:return None
        return offsets[group]+(i%2)*half+a%2
    E=[[(a,row(i,a))for a in range(4) if row(i,a)is not None]for i in range(6)]
    assert all(len({r for a,r in S})==len(S)for S in E)
    assert len({r for S in E for a,r in S})==sum(sizes)==sum(map(len,E))
    W={}
    for w in it.product(*E):
        a,r=tuple(zip(*w))
        if X(a)!=0 or any(a[i]//2 != a[i+1]//2 for i in (0,2,4)):continue
        assert len(set(r))==6 and r not in W
        count=sorted(sum(offsets[g]<=b<offsets[g]+sizes[g]for b in r)for g in range(6))
        assert count==[0,0,0,2,2,2]
        W[r]=a
    # Verify that all common positions are aligned and at a two-slot erasure
    # the symbol changes are precisely the common nonzero BCH phase delta.
    for p,q in it.combinations(W,2):
        bad=[i for i,(x,y) in enumerate(zip(p,q))if x!=y]
        assert len(set(p)&set(q))==6-len(bad) and len(bad)>=2
        if len(bad)==2:
            i,j=bad;a,b=W[p],W[q]
            assert a[i]^b[i]==a[j]^b[j]!=0
            assert X(phase(k,a[k]) for k in range(6))^X(phase(k,b[k])for k in range(6))==phase(i,a[i]^b[i])^phase(j,a[j]^b[j])
    # Old flat layer takes one block per group, so intersections are ≤3.
    inc=len(W)*(1<<13)
    assert len(W)==(256 if dim==96 else 160)
    records.append({'dimension':dim,'previous':baseline[dim],'new_supports':len(W),'increase':inc,'lower_bound':baseline[dim]+inc})
assert records==json.loads((P/'claims.json').read_text())
for r in records:print(r['dimension'],r['lower_bound'],'+',r['increase'])
