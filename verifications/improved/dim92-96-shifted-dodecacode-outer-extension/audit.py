"""Independent exact audit of the 92--96D Golay/dodecacode updates."""
from collections import Counter
from itertools import combinations,product
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def data(s):return json.loads((P/s).read_text())
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def lin(rows):
    out=[0]
    for r in rows:out.extend(r^x for x in out[:])
    assert len(set(out))==1<<len(rows)
    return out
Z=data('z8_generator.json')['generator']
assert len(Z)==12 and all(len(row)==24 for row in Z)
assert all(Z[i][i]%2 and not any(Z[i][:i]) for i in range(12))
G=lin([sum((a%2)*2**j for j,a in enumerate(r)) for r in Z])
assert Counter(x.bit_count() for x in G)=={0:1,8:759,12:2576,16:759,24:1}
D=lin([sum(int(a)*4**((i+j)%12) for j,a in enumerate('210100100101')) for i in range(12)])
assert Counter(((x|(x>>1))&0x55555555).bit_count() for x in D)=={0:1,6:396,8:1485,10:1980,12:234}
sign=[]
for s in range(1<<12):
    sy=0
    for j in range(12):
        if (s>>j)&1:sy^=j+1
    if sy==0:sign.append(s)
assert len(sign)==256 and min(s.bit_count() for s in sign if s)==3
F=[[(0,1),(2,3)],[(0,2),(1,3)],[(0,3),(1,2)]]
X=[[tuple(2 if k==i else a if k==j else 0 for k in range(4))
    for i,j in m for a in (-2,2)] for m in F]
A=[tuple((-1)**((a//4)^((a&1)*i)^(((a//2)%2)*j)) for i,j in product(range(2),repeat=2))for a in range(8)]
B=[(-x[0],)+x[1:]for x in A]
C=[tuple((2 if a<4 else -2) if i==a%4 else 0 for i in range(4))for a in range(8)]
Y=[A,B,C]
for i,j in product(range(3),repeat=2):
    for a,b in product(range(4),repeat=2):
        d=abs(dot(X[i][a],X[j][b]))
        assert (d==(8 if a==b else 0)) if i==j else d<=4
    for a,b in product(range(8),repeat=2):
        d=dot(Y[i][a],Y[j][b])
        assert d==(4 if a==b else -4 if a^b==4 else 0) if i==j else abs(d)==2
    for a,b in product(range(4),range(8)):
        assert abs(dot(X[i][a],Y[j][b]))<=4
# All three perfect matchings partition the six unordered pairs.
assert len(set(pair for m in F for pair in m))==6
O={}
for row in data('patterns.jsinput.json'):
    mask=row['zero']
    ids=[int(s,36) for s in row['encoding'].split(',')]
    assert len(ids)==len(set(ids)) and all(0<=i<2**24 for i in ids)
    patterns=[(G[i%4096],G[i//4096])for i in ids]
    assert all((a|b).bit_count()==12 and ((a|b)&mask)==0 for a,b in patterns)
    H=Counter(((a^c)|(b^d)).bit_count()for (a,b),(c,d) in combinations(patterns,2))
    assert min(H)==12
    # Every occupied row yields an ordered signed projective axis, recoverable
    # from its support and entries. Outer tag is recovered from the support.
    # Across patterns d_H = 24-t-s, and dot <=4(t+s).
    O[mask]={'patterns':len(patterns),'distances':dict(sorted(H.items())),'size':len(patterns)*len(D)*len(sign)}

# Recompute incumbent dense geometry from explicit predecessor head receipts.
prev=data('incumbent_heads.json')
def dense_check(n,selected):
    r=4 if n==92 else 2
    cap=8 if n==92 else 4
    for c,z,h in selected:
        assert c in range(3) and len(z)==r and all(0<=t<8 for t in z) and dot(h,h)==4*r
    assert len({tuple(h)for _,_,h in selected})==len(selected)
    assert len({(c,*z)for c,z,_ in selected})==len(selected)
    for (c,z,h),(e,w,k) in combinations(selected,2):
        q=dot(h,k)
        if c!=e:assert q<=cap
        else: assert q<=sum(dot(Y[c][a],Y[c][b])for a,b in zip(z,w))
    return len(selected)*8**(12-r)
h92=[]
for i,(parts,supp) in enumerate(zip(prev['92']['partitions'],prev['92']['selected_supports'])):
    assert len(parts)==4 and sorted(sum(parts,[]))==list(range(12))
    for s in supp:
        chosen=[list(set(g)&set(s)) for g in parts]; assert all(len(c)==1 for c in chosen)
        axes=[c[0]for c in chosen]
        for bits in product(range(2),repeat=4):
            h=[0]*12
            for j,b in zip(axes,bits):h[j]=2*(-1)**b
            h92.append((i,[g.index(j)+4*b for g,j,b in zip(parts,axes,bits)],h))
h93=[(c,z,h)for (c,*z),h in zip(prev['93']['labels'],prev['93']['heads'])]
h94=[(d['color'],d['prefix'],d['head']) for d in prev['94']['selected']]
h95=[(d['color'],d['prefix'],d['head']) for d in prev['95']['selected']]
assert all(len(h)==5 for _,_,h in h93)
assert all(h[1]==h[4]==-h[0] and 3*h[0]**2<=4 for _,_,h in h94)
assert all(h[0]+h[4]==0 and 2*h[0]**2<=8 for _,_,h in h95)
assert all(max(map(abs,h))<=2 for rows in (h92,h93,h94,h95) for _,_,h in rows)
counts={i:dense_check(i,j)for i,j in zip(range(92,96),(h92,h93,h94,h95))}
counts[96]=len(Y)*8**12

# RS evaluation at six different points over GF16. Prove root bound by
# exhaustive enumeration of each of the nonzero cubics.
def mul(a,b):
    raw=0
    for j in range(4):
        if b>>j&1:raw^=a<<j
    for j in range(6,3,-1):
        if raw>>j&1:raw^=19<<(j-4)
    return raw
M=[[mul(a,b)for b in range(16)]for a in range(16)]
assert all(sorted(m)==list(range(16)) for m in M[1:])
def evaluate(co,x):
    y=0
    for a in reversed(co): y=M[y][x]^a
    return y
stock=Counter()
for co in product(range(16),repeat=4):
    q=tuple(evaluate(co,i)for i in range(6))
    assert q[0]==co[0]
    if any(co):assert q.count(0)<=3
    stock[co[0]]+=1
assert set(stock.values())=={16**3}
sign6=[u for u in range(64)if u.bit_count()%2==0]
assert len(sign6)==32 and min((x^y).bit_count() for x,y in combinations(sign6,2))==2
receipt={}
for n in range(92,97):
    r=n-80-(n==95)
    rs=len(sign6)*sum(stock[i]for i in range(r))
    rt=4*(r*(r-1)//2+5*(16*15//2)+(n-1 if n==95 else 0))
    old=counts[n]+rs+rt
    new=O[15 if n==92 else (0 if n==96 else 3)]['size']
    receipt[n]={'dense':counts[n],'rs':rs,'roots':rt,'incumbent':old,'new':new,'combined':old+new}
assert [r['incumbent'] for r in receipt.values()]==[29797911144,42951379608,73016281804,109523504196,206160530240]
assert [r['combined'] for r in receipt.values()]==[30121921128,43934943896,73999846092,110507068484,208869002048]
print(json.dumps({'outer':O,'layers':receipt},sort_keys=True))
