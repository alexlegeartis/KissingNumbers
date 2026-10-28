"""Independent proof using shortened dual checks and packed support products.

Does not import replay.py. Z4 words are enumerated with six-plus-six halves;
small errors come from a binary Hamming ball. Sparse sections use expanded
integer signs and Reed--Solomon words use polynomial evaluation.
"""
from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
bits = getattr(int, 'bit_count', lambda x: bin(x).count('1'))
LUT16 = np.array([bits(i) for i in range(65536)], dtype=np.uint8)


def count_bits(a):
    a = np.asarray(a, dtype=np.uint64)
    if hasattr(np, 'bitwise_count'):
        return np.bitwise_count(a)
    return sum(LUT16[((a >> np.uint64(s)) & np.uint64(65535)).astype(np.uint16)] for s in (0, 16, 32, 48))


def binary_matrix(G):
    coefficients = ((np.arange(2**len(G), dtype=np.uint32)[:, None] >> np.arange(len(G))) & 1).astype(np.int16)
    B = coefficients @ G.astype(np.int16) % 2
    assert np.count_nonzero(B.sum(axis=1) == 0) == 1
    return B


def shortened_checks(G, inside):
    # Work on the complete 24-column generators. Eliminate the eight inside
    # coordinates, keeping four order-four rows and twice an order-four row.
    T = G.astype(np.int64).copy()
    at = 0
    for col in inside:
        candidates = [i for i in range(at, len(T)) if T[i, col] % 2]
        if not candidates:
            continue
        i = candidates[0]
        T[[at, i]] = T[[i, at]]
        T[at] = T[at]*int(T[at, col]) % 4
        for j in range(len(T)):
            if j != at:
                T[j] = (T[j]-T[j, col]*T[at]) % 4
        at += 1
    assert at == 7
    candidates = [(col, i) for col in inside for i in range(at, len(T)) if T[i, col]]
    col, i = candidates[0]
    T[[at, i]] = T[[i, at]]
    assert T[at, col] == 2
    for j in range(at+1, len(T)):
        if T[j, col]:
            T[j] = (T[j]-T[at]) % 4
    H = np.vstack([T[at+1:], 2*T[at]]) % 4
    assert H.shape == (5, 24) and not np.any(H[:, inside])
    assert not np.any(H @ G.T % 4)
    return H


def parent_proof(raw):
    G = np.array(raw['generator'], dtype=np.int64)
    assert G.shape == (12, 24) and np.all((G >= 0) & (G <= 7))
    assert not np.any(G @ G.T % 8)
    assert not np.any(np.tril(G[:, :12], -1)) and all(int(G[i, i]) % 2 for i in range(12))
    grey = np.array([[(a >> 2) & 1, ((a >> 2) ^ (a >> 1)) & 1,
                      ((a >> 2) ^ a) & 1, ((a >> 2) ^ (a >> 1) ^ a) & 1]
                     for a in range(8)], dtype=np.uint8)
    for a in range(8):
        assert sum(grey[a]) % 2 == 0
        for b in range(8):
            expected = 0 if a == b else 4 if (a-b) % 8 == 4 else 2
            assert sum(grey[a] != grey[b]) == expected
    G4 = G % 4
    digits = ((np.arange(4096)[:, None] >> (2*np.arange(6))) & 3).astype(np.int16)
    first = digits @ G4[:6] % 4
    second = digits @ G4[6:] % 4
    histogram = np.zeros(49, dtype=np.int64)
    for value in second:
        V = (first+value) % 4
        # Expanded binary Gray coordinates, not a Lee-weight lookup table.
        weight = np.sum(V >> 1, axis=1)+np.sum((V >> 1) ^ (V & 1), axis=1)
        histogram += np.bincount(weight, minlength=49)
    assert histogram.sum() == 2**24 and histogram[0] == 1 and not any(histogram[1:12])
    assert histogram[12] > 0
    packed = [sum(int(x % 2) << j for j, x in enumerate(row)) for row in G]
    residues = [0]
    for row in packed:
        residues += [word ^ row for word in residues]
    assert len(set(residues)) == 4096
    residue_histogram = [0]*25
    for word in residues:
        residue_histogram[bits(word)] += 1
    assert {bits(word) for word in residues} == {0, 8, 12, 16, 24}
    # Invert the ordinary Gray map on every binary error of weight <=3.
    ball = [0]+[sum(1 << c for c in support)
                for w in (1, 2, 3) for support in combinations(range(32), w)]
    inverse = np.array([0, 1, 3, 2], dtype=np.int64)
    errors = inverse[((np.array(ball, dtype=np.uint64)[:, None] >> (2*np.arange(16, dtype=np.uint64))) & 3).astype(np.uint8)]
    assert len(errors) == 5489
    octads = 0
    for label, word in enumerate(residues):
        if bits(word) != 8:
            continue
        octads += 1
        inside = [i for i in range(24) if word >> i & 1]
        outside = [i for i in range(24) if not (word >> i & 1)]
        H = shortened_checks(G4, inside)[:, outside]
        baseline = sum((G[i] for i in range(12) if label >> i & 1), np.zeros(24, dtype=np.int64)) % 8
        assert not np.any(baseline[outside] % 2)
        target = baseline[outside]//2
        wanted = H @ target % 4
        actual = errors @ H.T % 4
        assert np.all(np.any(actual != wanted, axis=1))
    assert octads == 759
    return grey, dict(z8_words=2**36, binary_union_words=2**37, binary_distance=24,
        z4_words_enumerated=2**24, z4_lee_histogram=histogram.tolist(),
        binary_residue_histogram=residue_histogram, octads_checked=759,
        low_lee_errors_per_octad=5489, low_lee_coset_exclusions=759*5489,
        generator_sha256=hashlib.sha256(G.astype(np.uint8).tobytes()).hexdigest(),
        grey_map=grey.tolist(), independent_uniform_ring_coordinates=list(range(12)))


def dual_signs(G):
    A = G.astype(np.uint8).copy()
    pivots = []
    at = 0
    for col in range(24):
        candidates = np.flatnonzero(A[at:, col])
        if not len(candidates):
            continue
        i = at+int(candidates[0])
        A[[at, i]] = A[[i, at]]
        for j in range(len(A)):
            if j != at and A[j, col]:
                A[j] ^= A[at]
        pivots.append(col)
        at += 1
    assert at == 14
    free = [i for i in range(24) if i not in pivots]
    D = np.zeros((len(free), 24), dtype=np.uint8)
    for i, col in enumerate(free):
        D[i, col] = 1
        D[i, pivots] = A[:, col]
    assert not np.any(D.astype(int) @ G.T.astype(int) % 2)
    return binary_matrix(D)


M4 = np.array([[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]], dtype=np.uint8)


def support_proof(data):
    lo, hi = data['supports_low'], data['supports_high']
    assert lo.dtype == hi.dtype == np.uint64 and lo.ndim == hi.ndim == 1 and lo.shape == hi.shape
    masks = [int(a) | (int(b) << 64) for a, b in zip(lo, hi)]
    assert len(set(masks)) == len(masks) and all(0 <= x < 2**96 and bits(x) == 24 for x in masks)
    maximum = pairs = 0
    for i in range(0, len(lo), 256):
        stop_i = min(i+256, len(lo))
        for j in range(i, len(lo), 256):
            stop_j = min(j+256, len(lo))
            products = count_bits(lo[i:stop_i, None] & lo[None, j:stop_j])+count_bits(hi[i:stop_i, None] & hi[None, j:stop_j])
            if i == j:
                products = products[np.triu_indices(stop_i-i, 1)]
            if products.size:
                assert np.all(products <= 12)
                maximum = max(maximum, int(products.max()))
                pairs += products.size
    assert pairs == len(masks)*(len(masks)-1)//2
    Gq = data['quaternary_generator']
    assert Gq.shape == (9, 24) and Gq.dtype == np.uint8 and np.all(Gq < 4)
    labels = np.arange(4**9, dtype=np.uint32)
    Q = np.zeros((len(labels), 24), dtype=np.uint8)
    for i in range(9):
        Q ^= M4[((labels >> (2*i)) & 3)[:, None], Gq[i][None, :]]
    qweights = np.count_nonzero(Q, axis=1)
    assert np.count_nonzero(qweights == 0) == 1 and min(qweights[1:]) == 12
    witness = data['dual_witnesses']
    assert witness.shape == (len(masks), 2) and witness.dtype == np.uint64 and np.all(witness < 2**48)
    occupied = [sum(bool(mask & (15 << (4*r))) for r in range(24)) for mask in masks]
    checks = 0
    for index, mask in enumerate(masks):
        used = set()
        misses = 0
        for encoded in witness[index]:
            encoded = int(encoded)
            if not encoded:
                continue
            d = [(encoded >> (2*r)) & 3 for r in range(24)]
            active = {r for r in range(24) if d[r]}
            assert not used & active
            used |= active
            for row in Gq:
                value = 0
                for a, b in zip(d, row):
                    value ^= int(M4[a, b])
                assert value == 0
            parity = 0
            for r in active:
                values = {int(M4[d[r], c]) >> 1 for c in range(4) if mask >> (4*r+c) & 1}
                assert len(values) == 1
                parity ^= values.pop()
            assert parity == 1
            misses += 1
        assert occupied[index]-misses <= 12
        checks += misses
    Gs = data['sign_generator']
    assert Gs.shape == (14, 24) and np.all(Gs <= 1)
    B = binary_matrix(Gs)
    assert min(B.sum(axis=1)[1:]) == 6
    D = dual_signs(Gs)
    assert min(D.sum(axis=1)[1:]) == 8
    assert np.any(D.sum(axis=1) == 8)
    retained = [0 if j % 2 else comb(j, j//2)*2**(14-j) for j in range(8)]+[70*2**7]
    histogram = [0]*9
    for mask in masks:
        histogram[bits(mask & 255)] += 1
    report = dict(supports=len(masks), support_pairs_checked=pairs,
        maximum_support_intersection=maximum, disjoint_dual_checks=checks,
        occupied_row_histogram=np.bincount(occupied, minlength=15).tolist(),
        sign_words=len(B), sign_distance=6, quaternary_words=len(Q), quaternary_distance=12,
        section95_support_points=sum(a*b for a,b in zip(histogram, retained)))
    return masks, B, report


def section_sign_counts(masks, B):
    signs = 1-2*B.astype(np.int16)
    prefix = np.column_stack([np.zeros(len(signs), dtype=np.int16), np.cumsum(signs, axis=1, dtype=np.int16)])
    profiles = Counter(tuple(bits((mask >> (4*r)) & 15) for r in range(8)) for mask in masks)
    zeros, equals = [0]*9, [0]*5
    for profile, count in profiles.items():
        offsets = [0]
        for n in profile:
            offsets.append(offsets[-1]+n)
        row_sums = np.column_stack([prefix[:, offsets[i+1]]-prefix[:, offsets[i]] for i in range(8)])
        for n in range(1, 9):
            zeros[n] += count*int(np.count_nonzero(np.all(row_sums[:, :n] == 0, axis=1)))
        for n in range(1, 5):
            equals[n] += count*int(np.count_nonzero(np.all(row_sums[:, :n] == row_sums[:, :1], axis=1)))
    return dict(occupancy_profiles=len(profiles), zero_counts=zeros, equal_counts=equals)


def field(q):
    def add(a, b):
        if q in (4, 16):
            return a ^ b
        if q == 9:
            return (a % 3+b % 3) % 3+3*((a//3+b//3) % 3)
        return (a+b) % q
    def mul(a, b):
        if q == 9:
            ar, ai, br, bi = a % 3, a//3, b % 3, b//3
            return (ar*br+2*ai*bi) % 3+3*((ar*bi+ai*br) % 3)
        if q in (4, 16):
            out = 0
            for j in range(4):
                if b >> j & 1:
                    out ^= a << j
            polynomial = 7 if q == 4 else 19
            for j in range(out.bit_length()-1, polynomial.bit_length()-2, -1):
                if out >> j & 1:
                    out ^= polynomial << (j-polynomial.bit_length()+1)
            return out
        return a*b % q
    A = np.array([[add(a,b) for b in range(q)] for a in range(q)], dtype=np.uint8)
    M = np.array([[mul(a,b) for b in range(q)] for a in range(q)], dtype=np.uint8)
    assert all(len(set(map(int, M[a]))) == q for a in range(1,q))
    for a,b,c in product(range(q), repeat=3):
        assert M[a, A[b,c]] == A[M[a,b], M[a,c]]
    return A,M


def small_layers():
    answer = {}
    for q in (9,11,13,16):
        A,M = field(q)
        labels = np.arange(q**4)
        coefficients = np.column_stack([(labels//q**j) % q for j in range(4)]).astype(np.uint8)
        W = np.zeros((len(labels),6), dtype=np.uint8)
        for j in range(3,-1,-1):
            W = A[M[W,np.arange(1,7)], coefficients[:,[j]]]
        weights = np.count_nonzero(W,axis=1)
        assert np.count_nonzero(weights==0)==1 and min(weights[1:])==3
        answer[str(q)] = dict(supports=q**4, signs=32, points=q**4*32, distance=3)
        if q == 16:
            selected = [2,2,1,1,1,1]
            hist = np.bincount(np.sum(W < np.array(selected), axis=1), minlength=7)
            S = np.array([[1-2*((word >> j)&1) for j in range(6)]
                          for word in range(64) if bits(word)%2],dtype=np.int16)
            assert len(S)==32 and np.all(np.diag(S@S.T)==6)
            dots=S@S.T;np.fill_diagonal(dots,-6);assert int(dots.max())==2
            retained=[int(np.count_nonzero(S[:,:j].sum(axis=1)==0)) for j in range(7)]
            answer['section95']=dict(selected_counts=selected,histogram=hist.tolist(),points=int(hist@np.array(retained)))
    return answer


def independent_bounds(grey, support, sparse, rs):
    even = Counter(map(int,(1-2*grey.astype(np.int16)).sum(axis=1)))
    odd_words = grey.copy();odd_words[:,0]^=1
    odd = Counter(map(int,(1-2*odd_words.astype(np.int16)).sum(axis=1)))
    result={}
    for dimension in range(88,97):
        constraints=[]
        if dimension<=92:
            m=96-dimension
            dense=even[0]**m*8**(12-m)
            transversal=0;extra=sparse['zero_counts'][m]
            q=9 if m==8 else 11 if m>=5 else 13
            small=rs[str(q)]['points'];pure=12*m+2*(96-4*m)
            text='each of the first m four-cell row sums is zero'
            for r in range(m):
                row=[0]*96;row[4*r:4*r+4]=[1]*4;constraints.append(row)
        elif dimension<=94:
            m=97-dimension
            dense=(sum(v**m for v in even.values())+sum(v**m for v in odd.values()))*8**(12-m)
            transversal=4**9*(2**14//2**(m-1));extra=sparse['equal_counts'][m]
            q=13;small=rs[str(q)]['points'];pure=12*m+2*(96-4*m)
            text='the first m four-cell row sums are equal'
            for r in range(1,m):
                row=[0]*96;row[:4]=[-1]*4;row[4*r:4*r+4]=[1]*4;constraints.append(row)
        elif dimension==95:
            m=2
            dense=(sum(v*even[-k] for k,v in even.items())+sum(v*odd[-k] for k,v in odd.items()))*8**10
            transversal=4**9*8192;extra=support['section95_support_points']
            small=rs['section95']['points'];pure=2*88
            text='the sum of the first eight coordinates is zero'
            constraints=[[1]*8+[0]*88]
        else:
            m=0;dense=2*8**12;transversal=4**9*16384;extra=support['supports']*16384
            small=16**4*32;pure=192;text='no section'
        # A nonzero minor modulo three certifies the real codimension.
        A=np.array(constraints,dtype=np.int64).reshape(96-dimension,96)%3
        rank=0
        for col in range(96):
            nz=np.flatnonzero(A[rank:,col])
            if not len(nz):continue
            k=rank+int(nz[0]);A[[rank,k]]=A[[k,rank]];A[rank]=A[rank]*int(A[rank,col])%3
            for i in range(len(A)):
                if i!=rank:A[i]=(A[i]-A[i,col]*A[rank])%3
            rank+=1
        assert rank==96-dimension
        if dimension<=94:
            assert 6*q<=96-4*m  # All weight-six points avoid the selected rows.
            roots=[]
            for a,b in combinations(range(4),2):
                for sign in (-1,1):
                    v=np.zeros(4,dtype=np.int16);v[a]=sign;v[b]=-sign;roots.append(v)
            R=np.array(roots);D=R@R.T;np.fill_diagonal(D,-2)
            assert len(R)==12 and np.all(R.sum(axis=1)==0) and int(D.max())==1
        assert 96-2*24<=48 and 4*(24-2*6)<=48 and 4*12<=48
        assert 4*2*6<=48 and 16*3<=48 and 16*(6-2*2)<=48
        assert 16*96<=48**2 and 16*48<=48**2
        terms=dict(dense=int(dense),transversal24=transversal,support24=extra,support6=small,pure=pure)
        result[str(dimension)]=dict(total=sum(terms.values()),terms=terms,section=text,
            selected_four_cell_rows=m,constraint_rows=constraints)
    return result


def main():
    if not __debug__:raise SystemExit('Run without -O: assertions are required.')
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    raw=json.loads((ROOT/'data/z8_generator.json').read_text())
    grey,parent=parent_proof(raw)
    with np.load(ROOT/'data/supports.npz',allow_pickle=False) as z:data=dict(z)
    masks,B,support=support_proof(data)
    sparse=section_sign_counts(masks,B);rs=small_layers()
    report=dict(status='PASS',representation='shortened dual checks and packed integer products',
        z8=parent,support=support,sparse_sections=sparse,reed_solomon=rs,
        dimensions=independent_bounds(grey,support,sparse,rs),
        support_sha256=hashlib.sha256((ROOT/'data/supports.npz').read_bytes()).hexdigest(),
        scope='Explicit compressed constructions and finite arithmetic proofs; no enumeration of 2^37 binary words. The parent code is published prior work.')
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',dimensions={d:r['total'] for d,r in report['dimensions'].items()})))


if __name__=='__main__':main()
