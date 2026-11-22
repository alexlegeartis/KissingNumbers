"""Replay the 46 full 93D bipartite-root Z8 fibres and retained body-only stock."""
from collections import Counter
from hashlib import sha1, sha256
from itertools import combinations, product
from pathlib import Path
import json

P=Path(__file__).resolve().parent
BASE="dim93-z8-variable-radius-heads-continuation"
FLAT="dim91-93-flat-parity-extensions"
TRIAL="dim91-93-triality-matching-transfer"
RS="dim89-96-rs-axis-code-repack"
STEINER="dim91-93-complementary-parity-steiner-append"
APPEND="dim91-93-nonsteiner-triality-append"
pinned={
 (BASE,"heads.json"):"a76ef2de11881e0edf83ba465b03b26b42cf575e",
 (BASE,"z8_generator.json"):"b9b7faca47e39c04091a73098f7e82dd5c76b065",
 (BASE,"verify.py"):"ca4233920c6324e6125e877b57404554d55a5f33",
 (BASE,"README.md"):"145701ae1ffa0fac52902f98cefcd1db1860ff37",
 (FLAT,"claims.json"):"d078130ffe70bfbc7383eca3c293c6ff9280a23f",
 (FLAT,"README.md"):"a5594b0abcca179c4219668df8013310d0f2ecf5",
 (TRIAL,"claims.json"):"ddfcafe618488c7e8cccdbd297ef356e020f436d",
 (TRIAL,"README.md"):"22c6311197fb3386ef293a88498a4250abd7c3cd",
 (RS,"claims.json"):"1b78a7a50a09fcbe7e032369d6770a4c8ce395da",
 (RS,"README.md"):"8c51adb336f99560ba70dc07e5c6b05da326cb80",
 (STEINER,"claims.json"):"4055f8326f454ca528215c0d47a1e2ef9c43200e",
 (STEINER,"README.md"):"25a02113ca1fc964d29c2f2dbb4ca672608777de",
 (APPEND,"claims.json"):"36a5c26324daf7ca6943af71622a410394ba8369",
 (APPEND,"README.md"):"468e0f2c8429fb45dd4b15ed45f20cff9ebba06b",
}
def pin(k):
    raw=(P.parent.joinpath(*k)).read_bytes()
    assert sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest()==pinned[k],k
    return raw
def claim(folder):
    J=json.loads(pin((folder,"claims.json")))
    if isinstance(J,dict): J=J["claims"]
    return next(x for x in J if x["dimension"]==93)

def dot(a,b): return sum(x*y for x,y in zip(a,b))
def alphabets():
    A=[tuple((-1)**(((q>>2)+(q&1)*a+((q>>1)&1)*b)%2)
             for a,b in product(range(2),repeat=2)) for q in range(8)]
    B=[(-a[0],*a[1:]) for a in A]
    C=[tuple((2 if q<4 else -2) if i==q%4 else 0 for i in range(4)) for q in range(8)]
    return [A,B,C]

def main():
    assert __debug__
    for k in pinned:pin(k)
    G0=pin((BASE,"z8_generator.json"))
    assert sha256(G0).hexdigest()=="903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73"
    G=json.loads(G0)["generator"]
    assert len(G)==12 and all(len(row)==24 and row[i]%2 and not any(row[:i]) for i,row in enumerate(G))
    assert len({(a*G[0][0]%8,(a*G[0][1]+b*G[1][1])%8) for a,b in product(range(8),repeat=2)})==64
    T=alphabets();W=(0,2,2,2,4,2,2,2)
    for A in T:
        assert len(set(A))==8
        assert all(dot(A[u],A[v])==4-2*W[(u-v)%8] for u,v in product(range(8),repeat=2))
    assert all(abs(dot(x,y))==2 for A,B in combinations(T,2) for x,y in product(A,B))
    old=json.loads(pin((BASE,"heads.json")))
    assert len(old)==len({(x["color"],*x["prefix"]) for x in old})==44
    raw=[tuple(map(int,l.split())) for l in (P/"witness.txt").read_text().splitlines()]
    assert all(len(r)==7 and 0<=r[0]<3 and 0<=r[1]<64 for r in raw)
    H=[(c,divmod(p,8),h) for c,p,*h in raw]
    assert len(H)==len({(c,p) for c,p,_ in H})==46
    expected=[]
    for i,j,s,t in product(range(2),range(2,5),range(2),range(2)):
        h=[0]*5;h[i]=2-4*s;h[j]=2-4*t
        expected.append((0,(i+4*s,j-2+4*t),tuple(h)))
    assert {(c,p,tuple(h)) for c,p,h in H if c==0}==set(expected)
    assert Counter(dot(h,h) for _,_,h in H)=={8:40,1:6}
    # The homogeneous minimum 24 of the pinned parent is inherited.
    for c,p,h in H:
        assert 2*(40+dot(h,h))<=88+dot(h,h)
    pairs=0
    for (c,p,h),(d,q,k) in combinations(H,2):
        b=44 if c!=d else 48-sum(dot(T[c][u],T[c][v]) for u,v in zip(p,q))
        x=b+dot(h,k)
        assert x<=0 or 4*x*x<=(88+dot(h,h))*(88+dot(k,k))
        pairs+=1
    assert pairs==1035
    # Body-only retained families, with conservative mode-independent caps.
    # Matching rows have two entries +/-sqrt(2); their per-row dot cap
    # is 2sqrt(2), with square 8, against all three dense alphabets.
    U=[tuple(s if t==i else w if t==j else 0 for t in range(4))
       for i,j in combinations(range(4),2) for s,w in product((-2,2),repeat=2)]
    assert len(U)==24 and max(abs(dot(x,y)) for x,y in product(U,sum(T,[])))==4
    for w in (10,11): assert 4*(4*w)**2<=8*w*88         # Golay
    for w in (5,6,7,8): assert 4*8*w*w<=4*w*88          # matching, incl. RS repack
    assert 4*22**2<=28*88                               # retained 7-row flats
    assert 4*16**2<=16*88                               # triality modes, including appends
    assert 4*4**2<=8*88                                 # single-row D4 roots
    # The recent sparse package chain retained the 44 dense fibres and
    # either appended body-only vectors or repacked body-only supports.
    f,t,r,s,a=(claim(k) for k in (FLAT,TRIAL,RS,STEINER,APPEND))
    assert f["lower_bound"]==t["previous"] and t["lower_bound"]==r["previous"]
    assert r["lower_bound"]==s["previous"] and s["lower_bound"]==a["previous"]
    previous=a["lower_bound"]
    assert previous==47398174480
    old_dense=44*8**10
    dense=46*8**10
    stock=previous-old_dense
    assert stock==153534224
    # Add all D93 roots with at least one head axis (head coords 88..92).
    # Dense coordinates have magnitude <=2, hence mixed dot <=4.
    # A head-body root tests at most one body coordinate of any sparse
    # ray. The listed pairs are (maximum coordinate square, norm square).
    roots=[(i,j,a,b) for i,j in combinations(range(93),2) if j>=88
           for a,b in product((-1,1),repeat=2)]
    assert len(set(roots))==len(roots)==1800
    assert sum(i>=88 for i,j,a,b in roots)==40
    D=[{i:a,j:b} for i,j,a,b in roots]
    assert all(sum(a*v.get(i,0) for i,a in u.items())<=1 for u,v in combinations(D,2))
    assert all(max(abs(x) for x in h)<=2 and 4*4**2<=2*(88+dot(h,h)) for _,_,h in H)
    assert max(abs(x) for w in sum(T,[]) for x in w)==2
    for m2,n2 in [(4,80),(4,88),(2,20),(2,24),(2,28),(2,32),
                  (1,28),(4,16),(1,2)]:
        assert 4*m2<=2*n2
    out={"dimension":93,"previous":previous,"fibres":len(H),"dense":dense,
         "retained_sparse":stock,"new_head_incident_roots":len(roots),
         "gain":dense-old_dense+len(roots),"lower_bound":dense+stock+len(roots)}
    assert out==json.loads((P/"claims.json").read_text())
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
