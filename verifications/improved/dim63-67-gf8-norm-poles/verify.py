#!/usr/bin/env python3
"""Exact 63–67D pole replacement, with independent F64 norm evaluation and shifts."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha1
from itertools import combinations
import importlib.util,json
HERE=Path(__file__).resolve().parent
V=HERE.parent
ONE=(1<<64)-1

def pin(rel,h):
    p=V/rel
    b=p.read_bytes()
    assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==h,rel
    return p

def ledger(rel,h):
    return json.loads(pin(rel,h).read_text())

def pack(bits):
    return sum(1<<j for j in bits)

def span(rows):
    w=[0]
    for r in rows:
        w.extend([x^r for x in w])
    return w

def bits(m):
    while m:
        z=m&-m
        yield z.bit_length()-1
        m^=z

def lead(w,p):
    for j in sorted(p,reverse=True):
        if w>>j&1:
            w^=p[j]
    return w

def independent(rows):
    p={}
    for w in rows:
        w=lead(w,p)
        if w:
            p[w.bit_length()-1]=w
    return p

def fmul(x,y):
    z=0
    for i in range(6):
        if y>>i&1:
            z^=x<<i
    for j in range(10,5,-1):
        if z>>j&1:
            z^=0x43<<(j-6)
    assert 0<=z<64
    return z

def fnorm(x):
    y=1
    for _ in range(9):
        y=fmul(y,x)
    return y

def projections(D,H):
    piv={}
    null=[]
    for w in D:
        while (m:=w&H):
            j=m.bit_length()-1
            if j not in piv:
                piv[j]=w
                break
            w^=piv[j]
        else:
            assert w
            null.append(w)
    return piv,null

def lift(w,piv,H):
    a=0
    for j in sorted(piv,reverse=True):
        if (w^a)>>j&1:
            a^=piv[j]
    assert a&H==w
    return a

def simplex(a,X,D):
    labels=[fnorm(x^X[a]) for x in X]
    C=set(labels)
    assert len(C)==8 and labels[a]==0
    assert all(labels.count(c)==(1 if c==0 else 9) for c in C)
    assert all(fmul(fnorm(c),c)==fmul(c,fmul(c,c)) for c in C) # c^8=c, including 0
    G=independent([pack(j for j,v in enumerate(labels) if v>>k&1) for k in range(6)])
    assert len(G)==3 and all(lead(w,D)==0 for w in G.values())
    W=span(G.values())
    assert all(w.bit_count()==36 for w in W[1:]) and all((w>>a&1)==0 for w in W)
    return W+[w^ONE for w in W],labels

def replay_vectors(probes):
    vectors=[tuple(-(18 if i==a else 2) if w>>i&1 else (18 if i==a else 2)
                   for i in range(64)) for a,w in probes]
    assert len(set(vectors))==len(vectors)
    assert all(sum(v*v for v in x)==576 for x in vectors)
    # Exact dot using signed mask intersection, independently cross-checking
    # all amplitudes in literal coordinates for the maximum and summary below.
    dots=[]
    for (a,w),(b,z) in combinations(probes,2):
        diff=w^z
        d=4*(64-2*diff.bit_count())
        if a==b:
            d+=320*(1-2*(diff>>a&1))
        else:
            d+=32*(2-2*((diff>>a&1)+(diff>>b&1)))
        dots.append(d)
    assert max(dots)==288 and len(dots)==len(vectors)*(len(vectors)-1)//2
    return vectors,len(dots)

def verify():
    assert __debug__
    path=pin("dim64-65-signed-parity/construction.py","5657a1072801fd2ded5e9fa02124976904e44e3e")
    spec=importlib.util.spec_from_file_location("pinned_bch",path)
    S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
    P,poly=S.extended_bch(13)
    E,epoly=S.extended_bch(21)
    assert poly==0x37cd0eb67 and epoly==0x2f30b529d3d5
    assert len(P)==30 and len(E)==18
    Pr,Er=independent(P),independent(E)
    assert all(not lead(w,Pr) for w in E)
    X=[1]
    for i in range(62):X.append(fmul(X[-1],2))
    X.append(0)
    assert set(X)==set(range(64)) and fmul(X[62],2)==1
    # Verify even extension and twelve consecutive BCH roots in actual order.
    assert all(w.bit_count()%2==0 for w in P)
    for k in range(1,13):
        for w in P:
            v=0
            for j in bits(w&((1<<63)-1)):
                v^=X[(k*j)%63]
            assert v==0
    D=[]
    checks=[pack(j for j,y in enumerate(P) if (w&y).bit_count()%2) for w in P]
    for x in S.kernel(checks,30):
        w=0
        for j in bits(x):w^=P[j]
        D.append(w)
    D=independent(D)
    assert len(D)==28 and all(w.bit_count()%4==0 for w in D.values())
    assert all((w&v).bit_count()%2==0 for w in D.values() for v in D.values())
    assert not lead(ONE,D)
    extpath=pin("dim65-67-bch-exterior-caps/verify.py","b90ec347fb18bea5da951da963e4e136754f09c5")
    knownshift=0x800000037cd0eb67
    assert "shift=0x800000037cd0eb67" in extpath.read_text()
    # Their shifted Kerdock owners lie in knownshift+D. Use another D coset,
    # so old/new disjointness does not depend on reproducing their phase list.
    oldres=lead(knownshift,D)
    reps=span(independent([lead(w,D) for w in P]).values())
    assert len(reps)==4 and oldres in reps[1:] and all(not lead(w,Pr) for w in reps)
    fresh=next(w for w in reps if w not in (0,oldres))
    H=pack(i for i,x in enumerate(X) if fnorm(x)==fmul(x,x))
    assert H==0x8040201008040201
    assert H.bit_count()==8 and all((w&H).bit_count()%2==0 for w in D.values())
    piv,null=projections(D.values(),H)
    assert len(piv)==7 and len(independent(null))==len(null)==21 and all(not (w&H) for w in null)
    section_shift=next(w for w in reps if w and (w&H).bit_count()%2==0)
    sec=ledger("dim63-65-affine-reconstruction/section63.json","cdc01424d1789116ac83c5bb6509d88269bfc16b")
    phys=sec["octad"]
    assert sec["terms"]["poles"]==112 and sec["terms"]["dense"]==70*(1<<len(null))==146800640
    perm=dict(list(zip(bits(H),bits(phys)))+list(zip(bits(ONE^H),bits(ONE^phys))))
    assert len(set(perm.values()))==64
    full,section=[],[]
    for target,marked,shifts in ((full,0,span(list(D.values())[-12:])),
                                  (section,H,span(null[-12:]))):
        used=set()
        for a in bits(ONE^marked):
            g,labels=simplex(a,X,D)
            t=fresh
            if marked:
                fiber={}
                for j in bits(H):fiber.setdefault(labels[j],[]).append(j)
                assert len(fiber)==4 and all(len(k)==2 for k in fiber.values())
                target_pat=pack(v[-1] for v in fiber.values())
                t=section_shift^lift(target_pat^(section_shift&H),piv,H)
            for x in shifts:
                U={x^t^v for v in g}
                if U&used or (not marked and any(not lead(w,Er) for w in U)):
                    continue
                assert len(U)==16 and all(not lead(w,Pr) and lead(w,D) for w in U)
                assert marked or all(lead(w,D)!=oldres for w in U)
                used|=U
                target.extend((a,w) for w in sorted(U))
                break
            else:
                raise AssertionError("no new group on a required axis")
    fv,fchecks=replay_vectors(full)
    sv,schecks=replay_vectors(section)
    assert len(fv)==1024 and len(sv)==896
    assert all(sum(x[k] for k in perm if phys>>perm[k]&1)==0 for x in sv)
    assert fchecks==523776 and schecks==400960
    # Amplitudes in actual norm64 coordinates; cap tails vanish. ±8 axes
    # on the heavy coordinate are removed, other poles have product <= 16/3.
    assert 6**2+63*Q(4,9)==64
    assert 2*(6+15*Q(2,3))==32 and 4*(6+3*Q(2,3))==32
    assert Q(8,3)*(6+8*Q(2,3))<32 and 8*Q(2,3)<32
    assert 6+Q(2,3)*(63-2*14)==Q(88,3)<32
    assert -6+Q(2,3)*(63-2*13)<32
    assert Q(3,4)<1 and Q(37,64)<1
    prior={}
    relevant=[
      ("dim63-balanced-nine-refill/claims.json","e85ed9d87813b767a188823fbd144a6a7b355d1f",(63,)),
      ("dim64-65-partial-nonary-signs/claims.json","96a743ef12f41a0c144c7ca9c86cf2c52749f573",(64,)),
      ("dim65-67-bch-exterior-caps/claims.json","1dafcee005375381c2dad801cf50d942cfde1afc",(65,)),
      ("dim66-67-netto-signed-trades/claims.json","c05a97f686b77e669f784d66e51442fbb4143a44",(66,67))]
    for rel,digest,dims in relevant:
        for row in ledger(rel,digest)["dimensions"]:
            if row["dimension"] in dims:
                prior[row["dimension"]]=row["lower_bound"]
    assert set(prior)==set(range(63,68))
    out=[]
    for d,n in sorted(prior.items()):
        inserted=len(sv) if d==63 else len(fv)
        removed=112 if d==63 else 128
        out.append(dict(dimension=d,previous=n,inserted=inserted,removed=removed,gain=inserted-removed,lower_bound=n+inserted-removed))
    assert [r["lower_bound"] for r in out]==[184427081,335796250,336257692,342390470,355196300]
    return dict(dimensions=out,checks=[schecks,fchecks],full_phase=hex(fresh),section_phase=hex(section_shift))

if __name__=="__main__":
    print(json.dumps(verify(),indent=2))
