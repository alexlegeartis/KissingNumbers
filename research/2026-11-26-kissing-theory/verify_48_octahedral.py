import hashlib,itertools,json,sys
from pathlib import Path
raw=Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(raw).hexdigest()=="5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162"
d=json.loads(raw); B=d["five_coordinate_supports"]; Q=d["four_coordinate_supports"]
O=((0,1,2,3,9,8),(0,3,5,7,8,9),(4,1,2,6,8,9),(4,5,6,7,9,8))
P=((1,2,2),(2,1,-2),(2,-2,1))
A,C,K,H=531,519,253,21
D=2*A; N=D*D
G=[]
for o in O:
    m=[]
    for u,s,p in itertools.product((-1,1),(-1,1),P):
        x=[0]*10
        for i,val in zip(o,(u*A,*(s*K*z for z in p),C,H)): x[i]=val
        assert sum(z*z for z in x)==N and max(map(abs,x))<=A
        assert sum(abs(x[j]) for j in o[:4])+2*abs(x[o[5]])<=2*D
        m.append(x)
    G.append(m)
F=sum(G,[])
assert len(F)==len({tuple(x) for x in F})==48
assert max(sum(a*b for a,b in zip(x,y)) for x,y in itertools.combinations(F,2))==N//2
paid=[]
for i,q in enumerate(Q):
    for s in range(16):
        if any(sum(x[k]*(-1 if s>>j&1 else 1) for j,k in enumerate(q))>D for x in F):
            paid.append((i,s))
assert paid==[(5,s) for s in range(4)]
assert max(sum(abs(x[k]) for k in b) for i,m in enumerate(G) for x in m for j,b in enumerate(B) if i!=j)==1058<D
assert max(len(set(a)&set(b)) for a,b in itertools.combinations(B,2))<=2
assert max(len(set(q)&set(b)) for q in Q for b in B)<=2
base=200864+20+224-len(paid)+len(F)
assert base==201152 and base+61==201213 and base+4*16==201216
print("PASS: 48 compatible donors, 4 skeleton payments; the Q4 argument in the note excludes rebuilding 61 rows for these donors.")

N=10089
F=[]
for o in O:
    for hot,sign,p in itertools.product((1,-1),(1,-1),P):
        r=[sign*24*z for z in p]
        if hot<0: r[0]=-r[0]
        x=[0]*10
        for i,val in zip(o,(hot*50,*r,49,2)): x[i]=val
        assert sum(z*z for z in x)==N and 4*max(map(abs,x))**2<=N
        assert (sum(abs(x[i]) for i in o[:4])+2*abs(x[o[5]]))**2<=4*N
        F.append(x)
assert len({tuple(x) for x in F})==48
assert all(2*sum(a*b for a,b in zip(x,y))<=N for x,y in itertools.combinations(F,2))
paid=[]
for i,q in enumerate(Q):
    for s in range(16):
        for x in F:
            dot=sum(x[k]*(-1 if s>>j&1 else 1) for j,k in enumerate(q))
            if dot>0 and dot*dot>N:
                paid.append((i,s));break
assert paid==[(5,s) for s in range(4)]
print("PASS: mirrored 48-donor site counts, exact skeleton payments and separation.")
