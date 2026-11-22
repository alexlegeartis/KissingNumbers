#!/usr/bin/env python3
"""Exact fixed-parent 96D partial triality fibres, standard library only."""
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
from hashlib import sha1, sha256

P = Path(__file__).resolve().parent
OLD = P.parent / 'dim96-affine-phase-projections'
MIX = P.parent / 'dim96-mixed-triality-six-rows'
HEX = P.parent / 'dim94-96-hexary-parity-packing'
pc = int.bit_count

def pin(p, digest):
    b = p.read_bytes()
    assert sha256(b).hexdigest() == digest, p
    return b.decode()

def pin_blob(p, digest):
    b = p.read_bytes()
    assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest() == digest, p
    return b.decode()

def mask_indices(m):
    return tuple(i for i in range(24) if m >> i & 1)

def word_weight(w):
    return pc((w | (w >> 1)) & 0x555)

def image(x, forms):
    return tuple(pc(x & a) % 2 for a in forms)

def code():
    out = set()
    for x, y in product(range(16), range(8)):
        a = image(x >> 1, (1,2,3,4,5,6))
        b = (15 if x & 1 else 0) ^ (8 if x & 4 else 0) ^ (1 if x & 8 else 0)
        c = image(y, (7,3,4,6,2,1))
        out.add(sum(((u ^ (x & 1)) + 2 * (v ^ (b >> i & 1))) << (2*i)
                    for i, (u,v) in enumerate(zip(a,c))))
    assert len(out) == 128 and min(word_weight(x^y) for x,y in combinations(out,2)) == 3
    return out

def tap(w, I):
    return sum((w >> (2*i) & 3) << (2*j) for j,i in enumerate(I))

@lru_cache(None)
def projected(W, I):
    return frozenset(tap(w, I) for w in W)

def unpack(a,e,o,W, allow_flat=False):
    assert all(0 <= s < 1 << 24 for s in (a,e,o))
    assert not (a&e or a&o or e&o)
    m = a | e | o
    assert pc(m) == 6 and (allow_flat or 1 <= pc(a) <= 3)
    assert W and len(W) == len(set(W)) and all(0<=w<4096 for w in W)
    return a,e,o,m,mask_indices(m),tuple(sorted(W))

def cross(X,Y):
    a,e,o,m,I,W = X
    b,f,p,n,J,V = Y
    q = pc(m & n)
    h = mask_indices((a&b)|(e&f)|(o&p))
    d = (q + len(h) - 5) // 2
    if d <= 0:
        return 0
    L = projected(W, tuple(I.index(i) for i in h))
    R = projected(V, tuple(J.index(i) for i in h))
    assert all(word_weight(x^y) >= d for x in L for y in R), (a,e,o,b,f,p)
    return 1

def ordinary():
    one = set()
    paired = set()
    for a,b,c,d,e,f in product(range(4), repeat=6):
        if a^b^c^d^e^f:
            continue
        v = (a,b,c,d,e,f)
        one.add(sum(1 << (4*j + x) for j,x in enumerate(v)))
        if a//2 == b//2 and c//2 == d//2 and e//2 == f//2:
            paired.add(sum(1 << (8*(j//2)+4*(x//2)+2*(j%2)+x%2)
                           for j,x in enumerate(v)))
    assert len(one) == 1024 and len(paired) == 256 and not one & paired
    return one | paired

def main():
    assert __debug__
    Q = code()
    signs = [x for x in range(64) if not pc(x)%2]
    assert len(signs) == 32 and min(pc(s^t) for s,t in combinations(signs,2)) == 2
    axis = [tuple(2*s*(i==j) for i in range(4)) for j in range(4) for s in (-1,1)]
    flat = list(product((-1,1), repeat=4))
    even = [x for x in flat if x.count(-1)%2 == 0]
    odd = [x for x in flat if x.count(-1)%2]
    frames = [axis,even,odd]
    gram = lambda x,y: sum(a*b for a,b in zip(x,y))
    assert all(len(F)==len(set(F))==8 and all(gram(x,x)==4 for x in F) for F in frames)
    assert all(gram(x,y) in (-4,0) for F in frames for x,y in combinations(F,2))
    assert all(abs(gram(x,y))==2 for F,G in combinations(frames,2) for x,y in product(F,G))
    shift = [tuple(a if k==i else b if k==j else 0 for k in range(4))
             for i,j in combinations(range(4),2) for a,b in product((-2,2),repeat=2)]
    assert all(abs(gram(x,y))<=4 for F in frames for x,y in product(F,shift))
    assert 6*4 == 24 and 4*24**2 == 24*96
    U = ordinary()

    prior = pin_blob(MIX/'witness.txt','a00a1bdd3f82f40ed930a8b894465a3feca824ee')
    prior += pin_blob(MIX/'witness_extra.txt','7c2e5db39c7780c573356242bd8f12ff3a82762d')
    R = [unpack(*map(int,l.split()),Q) for l in prior.splitlines()]
    for l in pin_blob(HEX/'best6_96.txt','bef34cfc76b36934686b7714d20299e53c1ea4ca').splitlines():
        m,t = map(int,l.split())
        assert t&m == t
        R.extend([unpack(0,t,m^t,Q,True),unpack(0,m^t,t,Q,True)])
    for l in pin(OLD/'witness.txt','b206025aaae7805529d93cc6b7f931cda8e894c974e056c6a82146e221a007bd').splitlines():
        a,e,o,p,t,_ = l.split()
        perm = tuple(map(int,p))
        assert sorted(perm) == list(range(6)) and 0 <= int(t) < 4096
        R.append(unpack(int(a),int(e),int(o),[tap(w^int(t),perm) for w in Q]))
    assert len(R) == 1170+438+103
    V = []
    for l in (P/'witness.txt').read_text().splitlines():
        a,e,o,n,*W = map(int,l.split())
        assert len(W) == n
        assert all(word_weight(x^y) >= 3 for x,y in combinations(W,2))
        V.append(unpack(a,e,o,W))
    assert len({v[:3] for v in V}) == len(V)
    assert not {v[:3] for v in V} & {r[:3] for r in R}
    for X in V:
        a,e,o,m,I,W = X
        gs = [15 << (4*j) for j in range(6)]
        # Mixed dots with amplitude-four RS, flat RS, grouped roots, and ordinary/paired BCH.
        assert sum(2 if a&g else int(bool(m&g)) for g in gs) * 4 <= 24
        assert sum(4 if (e|o)&g else 2*int(bool(m&g)) for g in gs) <= 12
        assert all(pc(a&g)<=1 for g in gs) and 9*48<=24**2
        assert all(2*(pc(m&u)+pc((e|o)&u))<=12 for u in U)
    tested_old = sum(cross(X,Y) for X,Y in product(V,R))
    tested_new = sum(cross(X,Y) for X,Y in combinations(V,2))
    base = json.loads(pin(OLD/'claims.json','c06fe1697229351b09fce88b07d5ff1b46acaa4358baab2c217cfe6f33dd8334'))
    assert base['dimension']==96 and base['lower_bound']==209292294976
    n = sum(len(X[-1]) for X in V)
    result = dict(dimension=96,previous=base['lower_bound'],fibres=len(V),axes=n,signs_per_axis=len(signs),
                  active_old_pairs=tested_old,active_new_pairs=tested_new,gain=n*len(signs),lower_bound=base['lower_bound']+n*len(signs))
    assert result == json.loads((P/"claims.json").read_text())
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
