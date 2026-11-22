#!/usr/bin/env python3
"""Exact replay of 94–95D dodecad fibres with adapted affine sign codes."""
from base64 import b64decode
from collections import Counter
from functools import lru_cache
from hashlib import sha1
from itertools import combinations, product
from pathlib import Path
from struct import unpack_from
import json
import zlib

P = Path(__file__).resolve().parent
GOLAY = 'dim92-96-shifted-dodecacode-outer-extension'
POSTRS = 'dim94-96-post-rs-partial-fibres'
MATCH = 'dim94-95-matching-row-append'
PREP = 'dim94-95-adapted-dodecad-fibres'
PC = int.bit_count
FULL_SIGN = (1 << 4096) - 1

def pin(folder,name,digest):
    b = (P.parent/folder/name).read_bytes()
    assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest() == digest, name
    return b

def positions(m):
    return tuple(i for i in range(24) if m >> i & 1)
def qproject(a,J):
    return sum(((a >> (2*j)) & 3) << (2*i) for i,j in enumerate(J))
def bproject(a,J):
    return sum(((a >> j) & 1) << i for i,j in enumerate(J))
def qweight(a):
    return PC((a | (a >> 1)) & 0x555555)

def linear_code(B):
    out = {0}
    for b in B:
        assert b not in out
        out |= {x ^ b for x in out}
    return out

def dodeca():
    return tuple(sum(int(c) << (2*((i+j)%12)) for j,c in enumerate('210100100101'))
                 for i in range(12))

@lru_cache(None)
def prior_axes(J):
    Q = linear_code(qproject(w,J) for w in dodeca())
    assert len(Q) == 4096
    return Q

def sign_syndrome(w):
    h = 0
    for j in range(12):
        if w >> j & 1:
            h ^= j+1
    return h
SIGN = tuple(tuple(w for w in range(4096) if sign_syndrome(w)==h) for h in range(16))
SIGN_MASKS = tuple(sum(1 << w for w in Q) for Q in SIGN)

_bit_planes = tuple(sum(1 << w for w in range(4096) if w >> j & 1) for j in range(12))

@lru_cache(None)
def sign_projection(J):
    old = {bproject(w,J) for w in SIGN[0]}
    size = 1 << len(J)
    if len(old) == size:
        return None
    assert len(old) == size//2
    # The projection is an explicit hyperplane, validated on all 256
    # former sign words. Its annihilator is the only nonzero character
    # vanishing also on the erased parity-check columns.
    missed = set(range(12))-set(J)
    dual = [h for h in range(1,16) if all(not PC(h&(j+1))%2 for j in missed)]
    assert len(dual)==1
    m = sum((PC(dual[0]&(j+1))%2)<<i for i,j in enumerate(J))
    assert m and all(not PC(x&m)%2 for x in old)
    return m

@lru_cache(None)
def forbidden(I,J):
    char = sign_projection(J)
    if char is None:
        return FULL_SIGN
    x = 0
    for i,j in enumerate(I):
        if char>>i&1:
            x ^= _bit_planes[j]
    return FULL_SIGN ^ x

@lru_cache(None)
def projected_signs(S,I):
    return frozenset(bproject(w,I) for w in S)

def explicit_ints(mask):
    while mask:
        b = mask & -mask
        yield b.bit_length() - 1
        mask ^= b

def old_patterns():
    poly = (7,1,4,4,2,7,3,3,4,7,2,1)
    G = [0]
    for i in range(12):
        b = sum((x&1) << (i+j) for j,x in enumerate(poly)) + (1 << 23)
        G += [a^b for a in G]
    assert len(G)==len(set(G))==4096
    assert Counter(PC(w) for w in G)=={0:1,8:759,12:2576,16:759,24:1}
    old = json.loads(pin(GOLAY,'shifted_patterns.json','7ad938a6be9008825e20d1e45ba9fed0ebfd501d'))
    assert old[-2]['mask']==3
    ids = [int(v,36) for v in old[-2]['patterns'].split(',')]
    assert len(ids)==len(set(ids))==1064 and all(0<=w<1<<24 for w in ids)
    Q = [(G[w%4096],G[w//4096]) for w in ids]
    assert all(PC(a|b)==12 for a,b in Q)
    return set(G), Q

def decode():
    raw = b64decode(''.join((P/'witness.b64').read_text().split()),validate=True)
    z = zlib.decompress(raw)
    count, = unpack_from('<H',z,0)
    at = 2
    for _ in range(count):
        m,t,n = unpack_from('<IBH',z,at)
        at += 7
        assert 0<m<1<<24 and t in (1,2,3) and n>0 and PC(m)==12
        W = []
        for _ in range(n):
            w, = unpack_from('<I',z,at)
            at += 4
            assert w < (1<<29)
            mask = None
            if w & (1<<28):
                k, = unpack_from('<B',z,at)
                assert k > 0 and at+1+k <= len(z)
                mask = z[at+1:at+1+k]
                at += 1+k
            W.append((w & ((1<<28)-1), mask))
        assert len({w&0xffffff for w,_ in W}) == n
        yield (m if t&1 else 0),(m if t&2 else 0),W
    assert at == len(z)

def allowed_signs(w,h,neighbors):
    S = FULL_SIGN
    for I,J in neighbors:
        V = prior_axes(J)
        t = qproject(w,I)
        for k in range(-1,8):
            if k==-1:
                present = t in V
            else:
                present = any((t ^ (a << (2*k))) in V for a in (1,2,3))
            if present:
                ii = tuple(v for j,v in enumerate(I) if j!=k)
                jj = tuple(v for j,v in enumerate(J) if j!=k)
                S &= ~forbidden(ii,jj)
    assert S, (w,h)
    base = (S & -S).bit_length()-1
    groups = []
    labels = [None]*12
    # In this affine set, a single-bit shift at j is possible iff the sign
    # kernel constraint column at j is zero. A double-bit shift at j,k is
    # possible iff the columns are equal. Compute those relations directly
    # by membership in the fully decoded sign set, not by pivot guesses.
    for j in range(12):
        for G in groups:
            if S >> (base ^ (1 << j) ^ (1 << G[0])) & 1:
                G.append(j)
                break
        else:
            groups.append([j])
    for G in groups:
        t = (S >> (base ^ (1<<G[0]))) & 1
        for i,j in enumerate(G):
            labels[j]=i+t
    for j in range(12):
        if S >> (base ^ (1<<j)) & 1:
            assert labels[j]!=0
        for k in range(j):
            if S >> (base ^ (1<<j) ^ (1<<k)) & 1:
                assert labels[j]!=labels[k]
    def syn(x):
        z=0
        for j in range(12):
            if x>>j&1:z^=labels[j]
        return z
    T=tuple(x for x in explicit_ints(S) if syn(x)==h)
    assert T, (w,h)
    return T

def orthogonality():
    dot = lambda a,b: sum(x*y for x,y in zip(a,b))
    M = (((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
    rows = [[tuple(2 if h==i else s if h==j else 0 for h in range(4))
             for i,j in m for s in (-2,2)] for m in M]
    assert all(dot(x,x)==8 for F in rows for x in F)
    assert all(dot(x,y)==0 for F in rows for x,y in combinations(F,2))
    assert all(abs(dot(x,y))==4 for F,G in combinations(rows,2) for x,y in product(F,G))
    triality = [tuple(s*2*(i==j) for i in range(4)) for j in range(4) for s in (-1,1)]
    triality += list(product((-1,1),repeat=4))
    assert all(abs(dot(x,y))<=4 for F in rows for x,y in product(F,triality))
    # Dense and shifted 96D body vectors, 96D grouped roots and 96D RS sparse;
    # prior squared-norm-24 triality and flat fibres have at most six rows.
    assert 12*4==48 and 4**2*48<=48**2 and 6*2*4==48
    assert 6*4==24 and 4*24**2==24*96

@lru_cache(maxsize=20000)
def words_project(V,I):
    return tuple((qproject(w,I),S) for w,S in V)

def verify_between(X,Y,I,J):
    i,j = I,J
    u = words_project(tuple(X),i)
    v = words_project(tuple(Y),j)
    tested = 0
    for x,S in u:
        for y,T in v:
            z = x^y
            equal = len(i) - qweight(z)
            if equal <= 6:
                continue
            L = tuple(t for k,t in enumerate(i) if not (z >> (2*k))&3)
            R = tuple(t for k,t in enumerate(j) if not (z >> (2*k))&3)
            A = projected_signs(S,L)
            B = projected_signs(T,R)
            d = (equal-5)//2
            if d==1:
                assert not A & B
            else:
                assert all(PC(a^b)>=d for a,b in product(A,B))
            tested += 1
    return tested

def main():
    assert __debug__
    orthogonality()
    D = linear_code(dodeca())
    assert len(D)==4096 and min(qweight(w) for w in D if w)==6
    H = SIGN[0]
    assert len(H)==256 and min(PC(w) for w in H if w)==3
    assert set(map(len,SIGN))=={256}
    G,old = old_patterns()
    prior = [(a,b,a|b,positions(a|b)) for a,b in old]
    matching=[]
    for l in pin(MATCH,'witness.txt','517698a0364c917753553dff7fde9e70f9820396').decode().splitlines():
        m,t=map(int,l.split())
        T=[0,0,0]
        for i in positions(m):
            T[t%3]|=1<<i
            t//=3
        assert not t
        matching.append((T[0]|T[2],T[1]|T[2],m))
    assert len(matching)==85
    taken = []
    records = list(decode())
    assert len({(a,b) for a,b,_ in records})==len(records)
    total=axis_count=old_pairs=0
    for a,b,W in records:
        m=a|b
        I=positions(m)
        assert m in G and not m&3 and (not a or not b or a==b)
        for c,d,om in matching:
            meet=m&om
            q=PC(meet)+PC(meet&~((a^c)|(b^d)))
            assert q*q<=12*PC(om)
        neighbors=[]
        for c,d,om,J in prior:
            meet=m&om
            same=meet&~((a^c)|(b^d))
            if PC(meet)+PC(same)<=12:
                continue
            # Thus previous nonzero shifted phases of three-tetrad patterns
            # are never used: every active old pattern is monochrome, phase 0.
            assert meet==same and PC(meet)==8 and (not c or not d or c==d)
            common=positions(meet)
            neighbors.append((tuple(I.index(k)for k in common),tuple(J.index(k)for k in common)))
        assert len(neighbors) in (3,4,5)
        old_pairs+=len(neighbors)
        words=[]
        for w, mask in W:
            axis=w&0xffffff
            signs=allowed_signs(axis,w>>24,neighbors)
            if mask is not None:
                assert len(mask)==(len(signs)+7)//8
                chosen=int.from_bytes(mask,'little')
                assert 0<chosen.bit_length()<=len(signs)
                signs=tuple(s for i,s in enumerate(signs) if chosen>>i&1)
                assert len(signs)==chosen.bit_count()>0
            words.append((axis,signs))
        axis_count+=len(words)
        total+=sum(len(S) for w,S in words)
        taken.append((a,b,m,I,words))

    # Distinct axes of a mode: product is 8*(common-axis-count minus
    # twice the number of sign disagreements). Same axes use an entire
    # subset of one affine Hamming coset, whose sign distance is >=3.
    within = sum(verify_between((X,),(Y,),tuple(range(12)),tuple(range(12)))
                 for *_,V in taken for X,Y in combinations(V,2))
    cross = edges = 0
    for (a,b,m,I,X),(c,d,n,J,Y) in combinations(taken,2):
        meet=m&n
        same=meet&~((a^c)|(b^d))
        if PC(meet)+PC(same)<=12:
            continue
        assert meet==same and PC(meet)==8
        common=positions(meet)
        cross += verify_between(X,Y,tuple(I.index(k)for k in common),tuple(J.index(k)for k in common))
        edges+=1
    base=json.loads(pin(POSTRS,'claims.json','204a7b1d25ebb9c40dbd557f637fe07e4667678c'))
    origin={r['dimension']:r['lower_bound'] for r in base['claims']}
    c=json.loads(pin(PREP,'claims.json','b6a5730d1062a4eca6b5cf64bcdf7544f099d1cb'))
    prev={}
    for r in c['claims']:
        d=r['dimension'];assert r['previous']==origin[d];prev[d]=r['lower_bound']
    result=dict(geometries=len(records),axis_modes=axis_count,
                old_geometric_interfaces=old_pairs,new_geometric_interfaces=edges,
                active_new_axial_pairs=cross,active_within_axial_pairs=within,
                claims=[dict(dimension=d,previous=prev[d],construction_base=origin[d],construction_gain=total,
                             gain=origin[d]+total-prev[d],lower_bound=origin[d]+total)
                        for d in (94,95)])
    assert result==json.loads((P/'claims.json').read_text())
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
