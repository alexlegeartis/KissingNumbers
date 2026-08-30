#!/usr/bin/env python3
"""The axis layers of dimensions 91, 92 and 93, from the RECORD configuration of R^k.

    python axis_record.py                re-derive and report against the shipped layers
    python axis_record.py --write        re-derive and write data/poles_rec_<k>.npz
    python axis_record.py --write 20     just that k

WHY THIS FILE EXISTS.  scripts/axis_rotate.py takes the poles to be a rotated copy of the
CAP DIRECTIONS, so the layer can never exceed |Z|.  In eight of the dimensions here |Z| is
smaller than tau(k) -- the caps are the minimal vectors of Lambda_k, and Lambda_k is not the
best kissing configuration of R^k -- so a rotated copy is capped below the ceiling however
good the rotation.  Nothing requires the poles to
be a copy of the caps: ANY 60-degree code at 30 degrees from them is a legal layer, and the
paper's Proposition on this construction asks nothing else of W.  At k = 19 and 20 the record
configuration is available in the very frame the cap file is written in; at k = 21 it is
built in R^21 and carried in by an orthogonal frame of the axis space.  Measured, the three
layers go from 10648, 17380 and 27682 to what this file reports, and the ceiling if the
rational rotation lost nothing would be tau(k) = 11948, 19448 and 29768.

WHAT THE SHIPPED CAP SET IS.  data/lam19_W.npy and data/lam20_W.npy are, after dividing by
two, integer vectors of norm 8 in exactly n = k coordinates, of two shapes: 4*C(n,2) vectors
(+-2,+-2,0,...,0), and 128 sign patterns on each octad of a shortened Golay code, always
with an EVEN number of minus signs.  That is Lambda_k, and it is also the frame Cohn-Li
(arXiv:2411.04916) work in.

THE RECORD, IN THAT FRAME.  Take the octad sign patterns ODD instead.  The count is the
same and the code is still 60 degrees -- but now a vector s = ((-1)^c_1, ..., (-1)^c_n) can
be adjoined for every c in the DUAL of the octad code, because <s,u> = +-8 would need s to
agree with an odd pattern on an octad, and s has an even number of minus signs there.  Two
such vectors are compatible when their Hamming distance is at least ceil(n/4), so the extra
term is a maximum independent set in a Cayley graph on the dual code.  That gives
10668 + 1280 = 11948 = tau(19) (Ho, arXiv:2603.10425) and 17400 + 2048 = 19448 = tau(20).

ONE INEQUALITY DOES ALL THE WORK, and it is not the one the search was built for.  A pole
needs only THIRTY degrees from the caps, not sixty, and at k = 19 and 20 -- where the source
shares the cap file's frame -- both new families clear that with NO ROTATION AT ALL: in the
halved frame an odd-sign octad vector meets a cap at |X| <= 6 against a limit of 6, and a
+-1 vector at |X| <= 8 against a limit of 10.  Only the 4*C(n,2) pair vectors fail, and they
fail completely -- they ARE cap directions, so all 684 of them sit at zero degrees.  The
unrotated layer is therefore already 9984 + 1280 = 11264 at k = 19 against the 10648
shipped, with denominator 1.  A rotation moves the pair vectors off their own caps and is
worth the rest.

THE POLE CONDITION IS ONE-SIDED.  Two poles are far enough apart when <a,a'> <= 1/2 -- an
obtuse pair is legal, and -a is a legal pole beside a.  The extras include complementary
codewords, which sit at 180 degrees, so testing |<a,a'>| would reject the record
configuration itself.  The cap condition stays two-sided, because the cap set contains -z
with z.

MIXED NORMS.  The poles no longer all have the same length: the base has norm 8 and the
extras norm n.  Every test here is written with the two norms explicit -- 4X^2 <= 3 |a|^2 m
D^2 -- rather than folded into a single constant.

WHAT IS HEURISTIC AND WHAT IS NOT.  The rotation is a search, and the independent set behind
the extras is a search.  Nothing rests on either: scripts/verify_poles.py rebuilds the whole
source from the shipped cap file and the shipped code words, checks that it is a 60-degree
code, and checks every kept pole against every cap in exact integer arithmetic.
"""
import json
import os
import sys
import time
from math import gcd, isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from axis_rotate import DENOMS, cayley, rationalise                # noqa: E402

PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
KS = (19, 20, 21)
THR = np.sqrt(3.0) / 2.0


# --------------------------------------------------------------------------- the source
def halved(k):
    """the shipped cap directions, divided by two: norm 8 in exactly k coordinates"""
    W = np.load(os.path.join(DATA, 'lam%d_W.npy' % k)).astype(np.int64)
    W = W[:, np.abs(W).sum(0) > 0]
    assert (W % 2 == 0).all(), 'k = %d is not written in the halvable frame' % k
    V = W // 2
    assert set(map(int, (V * V).sum(1))) == set([8]), 'the halved directions are not norm 8'
    assert V.shape[1] == k, 'the frame has %d coordinates, not %d' % (V.shape[1], k)
    return V


def census(V):
    """(octad supports in sorted order, the minus-sign parity they carry, the pair vectors)"""
    sup, pairs = {}, []
    for r in V:
        nz = np.nonzero(r)[0]
        if len(nz) == 8:
            sup.setdefault(tuple(nz), []).append(int((r < 0).sum() % 2))
        elif len(nz) == 2:
            pairs.append(r)
        else:
            raise AssertionError('a direction with %d nonzero coordinates' % len(nz))
    par = set()
    for s, ps in sup.items():
        assert len(ps) == 128, 'octad %s carries %d sign patterns, not 128' % (s, len(ps))
        par |= set(ps)
    assert len(par) == 1, 'the octad sign patterns are of mixed parity: %s' % sorted(par)
    return sorted(sup), par.pop(), np.array(pairs, dtype=np.int64)


def odd_base(oct_, pairs, n):
    """the pair vectors, then every ODD sign pattern on every octad, in a fixed order"""
    rows = [r for r in pairs]
    for O in oct_:
        for mask in range(256):
            if bin(mask).count('1') % 2 == 0:
                continue
            r = np.zeros(n, dtype=np.int64)
            for j, i in enumerate(O):
                r[i] = 1 - 2 * ((mask >> j) & 1)
            rows.append(r)
    Z = np.array(rows, dtype=np.int64)
    assert set(map(int, (Z * Z).sum(1))) == set([8])
    return Z


def dual_code(oct_, n):
    """the dual of the octad code: every c whose intersection with every octad is even"""
    M = [sum(1 << (n - 1 - i) for i in O) for O in oct_]
    rows, piv = list(M), []
    for col in range(n - 1, -1, -1):
        b = 1 << col
        idx = next((i for i, a in enumerate(rows) if a & b), None)
        if idx is None:
            continue
        pr = rows.pop(idx)
        rows = [a ^ pr if a & b else a for a in rows]
        piv.append((col, pr))
    assert all(a == 0 for a in rows)
    pivcols = set(c for c, _ in piv)
    ker = []
    for f in [c for c in range(n) if c not in pivcols]:
        x = 1 << f
        for col, a in reversed(piv):
            if bin((a ^ (1 << col)) & x).count('1') % 2:
                x |= 1 << col
        for a in M:
            assert bin(a & x).count('1') % 2 == 0, 'a kernel word fails a parity check'
        ker.append(x)
    D = [0]
    for b in ker:
        D = D + [x ^ b for x in D]
    assert len(D) == 2 ** len(ker)
    return D, len(piv), ker


def vectors_of(code, n):
    """the +-1 vectors of a list of binary words"""
    return np.array([[1 - 2 * ((c >> (n - 1 - i)) & 1) for i in range(n)]
                     for c in sorted(code)], dtype=np.int64)


def caps(k):
    """the cap directions as the shipped file writes them, zero columns dropped"""
    W = np.load(os.path.join(DATA, 'lam%d_W.npy' % k)).astype(np.int64)
    return W[:, np.abs(W).sum(0) > 0]


def source(k, code=None, M=None):
    """(Zk, T, split, code, M, mu) -- the record configuration and the pulled-back caps

    Where the cap file is already written in k coordinates the octads are READ from it, so
    the odd base sits in the same frame as the caps and inherits the 30-degree condition for
    free.  Where it is not, the octads come from the Golay code and the configuration is
    built in R^k, then carried into the axis space by the frame M.
    """
    W = caps(k)
    if M is None:
        M, mu = frame(k, W)
        assert M is not None, 'k = %d: no orthogonal frame of the axis space' % k
    mu = int((M.T @ M)[0, 0])
    assert (M.T @ M == mu * np.eye(M.shape[1], dtype=np.int64)).all(), 'M is not a similarity'
    if W.shape[1] == M.shape[1]:
        V = halved(k)
        n = V.shape[1]
        oct_, parity, pairs = census(V)
    else:
        n = M.shape[1]
        oct_ = octads_shortened(n)
        pairs = pair_vectors(n)
        assert len(pairs) == 4 * (n * (n - 1) // 2)
    Z = odd_base(oct_, pairs, n)
    if code is None:
        code = extras_code(k, oct_, n)
    E = vectors_of(code, n)
    T = np.concatenate([Z, E]) if len(E) else Z
    Zk = W @ M                                   # the caps pulled back to R^k
    assert (np.einsum('ij,ij->i', Zk, Zk) == mu * int(W[0] @ W[0])).all(), \
        'the pull-back changed the cap norms: the frame does not span the axis space'
    return Zk, T, (len(pairs), len(Z), len(T)), sorted(code), M, mu, oct_, n



# ------------------------------------------------- building the source where it is not read
GOLAY_BASIS = os.path.join(HERE, '..', '..', '..', '..', 'common', 'data', 'golay_basis.txt')


def golay():
    """the 4096 codewords of the extended binary Golay code, as 24-bit integers"""
    B = [int(''.join(ln.split()), 2) for ln in open(GOLAY_BASIS) if ln.strip()]
    assert len(B) == 12, 'the Golay basis file has %d rows' % len(B)
    words = [0]
    for b in B:
        words = words + [w ^ b for w in words]
    assert len(set(words)) == 4096
    wt = sorted(bin(w).count('1') for w in words)
    assert wt[0] == 0 and wt.count(8) == 759 and set(wt) == set([0, 8, 12, 16, 24]), \
        'this is not the extended Golay code'
    return words


def octads_shortened(n):
    """the weight-8 Golay words supported inside the first n coordinates, as index tuples"""
    out = []
    for w in golay():
        if bin(w).count('1') != 8:
            continue
        idx = [i for i in range(24) if w >> (23 - i) & 1]
        if max(idx) < n:
            out.append(tuple(idx))
    return sorted(out)


def pair_vectors(n):
    """the 4 C(n,2) vectors (+-2, +-2, 0, ..., 0)"""
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            for a in (2, -2):
                for b in (2, -2):
                    r = np.zeros(n, dtype=np.int64)
                    r[i], r[j] = a, b
                    rows.append(r)
    return np.array(rows, dtype=np.int64)


def norm4_in_span(W):
    """the integer vectors of norm 4 of the ambient lattice that lie in the span of W"""
    A = W.astype(np.float64)
    _u, sv, vt = np.linalg.svd(A, full_matrices=False)
    r = int((sv > 1e-8 * sv[0]).sum())
    B = vt[:r]
    n = W.shape[1]
    rows = []
    for i in range(n):
        for sgn in (2, -2):
            v = np.zeros(n, dtype=np.int64)
            v[i] = sgn
            rows.append(v)
    import itertools
    for c in itertools.combinations(range(n), 4):
        for sg in itertools.product((1, -1), repeat=4):
            v = np.zeros(n, dtype=np.int64)
            for j, i in enumerate(c):
                v[i] = sg[j]
            rows.append(v)
    V = np.array(rows, dtype=np.int64)
    R = V.astype(np.float64) - (V.astype(np.float64) @ B.T) @ B
    return V[np.abs(R).max(1) < 1e-7], r


def frame(k, W, seed=0):
    """(M, mu): k columns spanning the axis space with M^T M = mu I, or (None, None)

    The columns are looked for among the norm-4 vectors of the ambient lattice that lie in
    the span, which is the class the discriminant allows at k = 21.  Where the span is the
    whole of R^k the frame is the identity and nothing is searched.
    """
    r = int(np.linalg.matrix_rank(W.astype(np.float64)))
    if r == W.shape[1]:
        return np.eye(r, dtype=np.int64), 1
    V4, r = norm4_in_span(W)
    rng = np.random.default_rng(seed)
    for _ in range(40000):
        chosen = []
        cand = np.arange(len(V4))
        while len(cand) and len(chosen) < r:
            i = int(cand[rng.integers(len(cand))])
            chosen.append(i)
            cand = cand[(V4[cand] @ V4[i]) == 0]
        if len(chosen) == r:
            M = V4[chosen].T
            G = M.T @ M
            assert (G == 4 * np.eye(r, dtype=np.int64)).all()
            return M, 4
    return None, None


# ------------------------------------------------------------------ the extra pole vectors
def nclasses(S0, H):
    return len(set(min(s ^ x for x in H) for s in S0))


def independent_sets(adj, m, lo):
    """every independent set of the quotient graph of size at least lo, as bitmasks"""
    out = []

    def go(i, cur, size, forbidden):
        if size + (m - i) < lo:
            return
        if i == m:
            if size >= lo:
                out.append(cur)
            return
        if not (forbidden >> i & 1):
            go(i + 1, cur | (1 << i), size + 1, forbidden | adj[i])
        go(i + 1, cur, size, forbidden)

    go(0, 0, 0, 0)
    return out


def quotient(K, H, S0):
    """(class representatives, adjacency bitmasks) of Cay(K, S0) modulo H"""
    S0set = set(S0)
    keys = sorted(set(min(v ^ x for x in H) for v in K))
    m = len(keys)
    adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if any((keys[i] ^ keys[j] ^ x) in S0set for x in H):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return keys, adj


def alpha_of(adj, m):
    """the independence number of the quotient, by branch and bound on bitmasks"""
    comp = [((1 << m) - 1) ^ adj[i] ^ (1 << i) for i in range(m)]
    best = [0]

    def go(size, cand):
        if size + bin(cand).count('1') <= best[0]:
            return
        if cand == 0:
            best[0] = max(best[0], size)
            return
        i = (cand & -cand).bit_length() - 1
        go(size + 1, cand & comp[i] & ~((1 << (i + 1)) - 1))
        go(size, cand ^ (1 << i))

    go(0, (1 << m) - 1)
    return best[0]


def best_pair(adj, m):
    """the largest |I1| + |I2| over DISJOINT independent sets of the quotient

    Listing every independent set and picking the big ones is not affordable: at the first
    beam depth the subgroup has order two and the quotient has 256 classes, where that list
    is astronomical.  Take the independence number first -- branch and bound, cheap -- and
    enumerate only the sets within one of it.
    """
    top = alpha_of(adj, m)
    for lo in (top, top - 1):
        if lo < 1:
            continue
        big = independent_sets(adj, m, lo)
        best = None
        for a in big:
            for b in big:
                if a & b == 0 and (best is None
                                   or bin(a).count('1') + bin(b).count('1') > best[0]):
                    best = (bin(a).count('1') + bin(b).count('1'), a, b)
        if best:
            return best
    return None


def split_beam(comp, S, width=3000, evaluate=200, verbose=True):
    """two disjoint coset unions across an index-two split; (size, the set)

    The connection set S has one word of weight below the rest.  Where that word w lies
    outside the subgroup K generated by the others, the component splits as K u (w + K), and
    a difference ACROSS the halves is in S only when it is w itself -- so the two halves need
    only be disjoint, not mutually independent.  Each half may then use cosets of the same
    subgroup H.

    CHOOSING H IS THE WHOLE GAME, and it is the choice, not the coset structure, that decides
    whether the record is reached.  What multiplies |H| is the independence number of the
    quotient, which is large only when the image of S is small, so H must absorb the pairwise
    differences of the connection words.  Taking the single best merge at each step stalls at
    eight classes at k = 19 -- 4 * 64 = 256 per component, exactly Cohn-Li's count and 256
    poles short of Ho's.  A beam over the same moves reaches five classes, the quotient is
    then the Clebsch graph on sixteen vertices, its independence number is 5, two disjoint
    such sets exist, and 2 * 5 * 32 = 320 per component is Ho's record.

    coset_beam() below reaches 320 as well, WITHOUT the split, with |H| = 64 and a quotient
    of independence 5.  So the split is not what buys the record; the guided choice of H is.
    Both routes are run and the larger is taken.
    """
    best = None
    for w in S:
        K = set([0])
        for s in S:
            if s != w and s not in K:
                K |= set(x ^ s for x in K)
        if w in K or len(K) * 2 != len(comp):
            continue
        S0 = [s for s in S if s != w]
        if not set(S0) <= K:
            continue
        S0set = set(S0)
        pool = sorted((set(a ^ b for a in S0 for b in S0) & K) - S0set - set([0]))
        beam, seen = [frozenset([0])], set()
        for depth in range(1, 7):
            nxt = []
            for H in beam:
                for g in pool:
                    if g in H:
                        continue
                    H2 = frozenset(H | set(x ^ g for x in H))
                    if H2 in seen or (H2 & S0set):
                        continue
                    seen.add(H2)
                    nxt.append((nclasses(S0, H2), H2))
            if not nxt:
                break
            nxt.sort(key=lambda t: t[0])
            beam = [h for _c, h in nxt[:width]]
            for H in beam[:evaluate]:
                keys, adj = quotient(K, H, S0)
                if len(keys) > 24:      # a quotient this large cannot be solved exactly, and
                    continue            # a subgroup this small cannot pay for it either
                got = best_pair(adj, len(keys))
                if not got:
                    continue
                size = got[0] * len(H)
                if best is None or size > best[0]:
                    A = [keys[i] ^ x for i in range(len(keys)) if got[1] >> i & 1 for x in H]
                    B = [w ^ keys[i] ^ x for i in range(len(keys)) if got[2] >> i & 1
                         for x in H]
                    best = (size, A + B, len(H), len(keys), nclasses(S0, H), got[0])
                    if verbose:
                        print('      |H| = %-4d quotient %2d classes, image of S is %d, '
                              'two disjoint sets of %d -> %d per component'
                              % (best[2], best[3], best[4], best[5] // 2, size))
    return best


def coset_beam(comp, S, width=3000, evaluate=200, verbose=True):
    """largest |H| * alpha(comp/H) over subgroups H of comp missing S; (size, the set)"""
    Sset = set(S)
    # The generators are drawn from the XOR-differences of S.  Adding the rest of the group
    # multiplies the beam expansion by a factor of five at k = 19 and 21 and finds nothing
    # the differences do not: what the beam is steering by is how far the image of S has
    # collapsed, and only a difference collapses it.
    pool = sorted((set(a ^ b for a in S for b in S) & set(comp)) - Sset - set([0]))
    if not pool:
        pool = [c for c in comp if c and c not in Sset]
    beam, seen = [frozenset([0])], set()
    best = None
    for depth in range(1, len(bin(len(comp))) - 2):
        nxt = []
        for H in beam:
            for g in pool:
                if g in H:
                    continue
                H2 = frozenset(H | set(x ^ g for x in H))
                if H2 in seen or (H2 & Sset):
                    continue
                seen.add(H2)
                nxt.append((nclasses(S, H2), H2))
        if not nxt:
            break
        nxt.sort(key=lambda t: t[0])
        beam = [h for _c, h in nxt[:width]]
        for H in beam[:evaluate]:
            keys, adj = quotient(comp, H, S)
            if len(keys) > 24:
                continue
            a = alpha_of(adj, len(keys))
            size = a * len(H)
            if best is None or size > best[0]:
                sel = independent_sets(adj, len(keys), a)[0]
                out = [keys[i] ^ x for i in range(len(keys)) if sel >> i & 1 for x in H]
                best = (size, out, len(H), len(keys), a)
                if verbose:
                    print('      |H| = %-4d quotient %2d classes, independence %d -> %d '
                          'per component' % (best[2], best[3], a, size))
    return best


def extras_code(k, oct_, n):
    """the binary words whose +-1 vectors may be adjoined to the odd base"""
    D, rank, ker = dual_code(oct_, n)
    dmin = (n + 3) // 4
    S = [d for d in D if d and bin(d).count('1') < dmin]
    comp = set([0])
    for s in S:
        if s not in comp:
            comp |= set(x ^ s for x in comp)
    reps = sorted(set(min(v ^ x for x in comp) for v in D))
    print('   k = %d: octad code dim %d, dual dim %d (%d words), minimum distance %d, '
          '|S| = %d' % (k, rank, len(ker), len(D), dmin, len(S)))
    print('      the graph is %d component(s) of %d' % (len(reps), len(comp)))
    got = split_beam(sorted(comp), S)
    direct = coset_beam(sorted(comp), S)
    if direct and (got is None or direct[0] > got[0]):
        got = direct
    assert got, 'neither route produced an independent set'
    out = sorted(set(r ^ x for r in reps for x in got[1]))
    assert len(out) == len(reps) * got[0]
    for i in range(len(out)):                       # the search is heuristic; this is not
        for j in range(i + 1, len(out)):
            assert bin(out[i] ^ out[j]).count('1') >= dmin, 'the code is not distance %d' % dmin
    print('      -> %d extra pole vectors' % len(out))
    return out


# ------------------------------------------------------------------------ the exact test
def cleared(T, V, N, D):
    """boolean per pole of T: does N/D carry it clear of every 30-degree cap?

    Exact throughout.  With X = (N t)^T z an integer, the condition 4<Qt,z>^2 <= 3|t|^2|z|^2
    is 4 X^2 <= 3 |t|^2 m D^2, and 3 |t|^2 m D^2 is never a perfect square here, so it is
    2|X| <= isqrt(3 |t|^2 m D^2) with one big-integer square root per norm and int64 after.
    """
    D = int(D)
    m = int(V[0] @ V[0])
    Ni = np.array(N.tolist(), dtype=np.int64)
    assert (Ni.astype(object) == N).all(), 'N does not fit in int64'
    A = T @ Ni.T                                    # D * (rotated poles)
    CZ = V.T
    # the two bounds scripts/verify_poles.py asserts while splitting A at bit 31
    _Ahi = A >> 31
    _bd = (int(np.abs(_Ahi).max()) + 2 ** 31) * int(np.abs(CZ).max()) * T.shape[1]
    assert _bd < 2 ** 50, 'D too large: verify_poles.py could not split this in float64'
    assert int(np.abs(_Ahi.astype(np.float64) @ CZ.astype(np.float64)).max()) < 2 ** 31, \
        'D too large: verify_poles.py could not shift the high limb'
    norms = (T * T).sum(1)
    lim = {}
    for na in sorted(set(int(v) for v in norms)):
        L = 3 * na * m * D * D
        lim[na] = isqrt(L)
        assert lim[na] ** 2 < L < (lim[na] + 1) ** 2, '3|t|^2 m D^2 is a perfect square'
    # TWO FLOAT LIMBS, not one int64 product.  numpy has no BLAS for int64, and this
    # product is 29768 x 21 by 21 x 27720 at k = 21 -- a minute and a half, which the search
    # pays once per candidate denominator.  Split A at bit 31 as verify_poles.py does: the
    # assert above pins every limb product under 2^50, so float64 holds it and every partial
    # sum exactly, and the recombination is int64.  Same numbers, fifty times faster.
    SH = 31
    Hf = (A >> SH).astype(np.float64)
    Lf = (A - ((A >> SH) << SH)).astype(np.float64)
    CZf = CZ.astype(np.float64)
    good = np.empty(len(T), bool)
    worst = 0
    B = max(1, int(8e6 // max(1, len(V))))
    for s in range(0, len(T), B):
        XH = Hf[s:s + B] @ CZf
        XL = Lf[s:s + B] @ CZf
        assert (XH == np.rint(XH)).all() and (XL == np.rint(XL)).all(), 'not integral'
        XHi = XH.astype(np.int64)
        assert int(np.abs(XHi).max()) < 2 ** 31, 'the shift would overflow'
        X = (XHi << SH) + XL.astype(np.int64)
        mx = np.abs(X).max(1)
        worst = max(worst, int(mx.max()))
        lo = np.array([lim[int(v)] for v in norms[s:s + B]], dtype=np.int64)
        good[s:s + B] = 2 * mx <= lo
    assert 2 * worst < 2 ** 63 - 1, 'the int64 product could overflow'
    return good


def is_code(T, split, report=True):
    """the 60-degree condition on the whole source, exactly and one-sided"""
    norms = (T * T).sum(1)
    worst, bad = -10 ** 9, 0
    B = max(1, int(4e6 // len(T)))
    for s in range(0, len(T), B):
        X = T[s:s + B] @ T.T
        for i in range(X.shape[0]):
            X[i, s + i] = -10 ** 9                  # a pole against itself is not a pair
        na = norms[s:s + B][:, None]
        worst = max(worst, int(X.max()))
        bad += int(((X > 0) & (4 * X * X > na * norms[None, :])).sum())
    if report:
        print('      the source is a 60-degree code: largest positive inner product %d, '
              '%d pairs break it' % (worst, bad))
    return bad == 0


# ------------------------------------------------------------------------- the rotation
def unit(A):
    return A / np.linalg.norm(A, axis=1, keepdims=True)


def active(Tn, Zn, Q, margin, cap=150000, blk=2000):
    """the (pole, cap) pairs within `margin` of the threshold"""
    TQ = Tn @ Q.T
    rows, cols = [], []
    for t in range(0, len(TQ), blk):
        r, c = np.nonzero(np.abs(TQ[t:t + blk] @ Zn.T) > THR - margin)
        rows.append(r + t)
        cols.append(c)
    r, c = np.concatenate(rows), np.concatenate(cols)
    if len(r) > cap:
        v = np.abs(np.einsum('ij,ij->i', TQ[r], Zn[c]))
        keep = np.argsort(-v)[:cap]
        r, c = r[keep], c[keep]
    return r, c


def minimax(Tn, Zn, Q, steps=900, eta=0.03, beta=200.0, margin=0.15, refresh=50):
    """push the worst cosine down, on an ACTIVE SET of near-threshold pairs

    The all-pairs objective that settled k = 7, 9 and 10 is 19448 x 17400 here -- 2.7 GFLOP a
    step -- but almost every pair is far below the threshold and contributes nothing to a
    soft argmax.  Keeping the pairs within a margin of it, and refreshing that set against
    all pairs every few dozen steps, is the same descent at a fiftieth of the cost.
    """
    r, c = active(Tn, Zn, Q, margin)
    best, bQ = float('inf'), Q.copy()
    for it in range(steps):
        if it % refresh == 0 and it:
            r, c = active(Tn, Zn, Q, margin)
            if len(r) == 0:
                break
        A, B = Tn[r] @ Q.T, Zn[c]
        s = np.einsum('ij,ij->i', A, B)
        f = np.abs(s)
        if f.max() < best:
            best, bQ = float(f.max()), Q.copy()
        w = np.exp(beta * (f - f.max()) / THR)
        w /= w.sum()
        M = (((w * np.sign(s))[:, None] * B).T @ Tn[r]) @ Q.T
        M = 0.5 * (M - M.T)
        ev, U = np.linalg.eigh(1j * (-eta * M))     # expm of a real skew matrix
        Q = np.real((U * np.exp(-1j * ev)) @ U.conj().T) @ Q
        u, _s, vt = np.linalg.svd(Q)
        Q = u @ vt
        if it % 300 == 299:
            eta *= 0.5
            beta *= 2.0
    return bQ


def refine(Q, Nb, Db, C, exact, Tf, Zf, m, denoms=None, keep=3):
    """(count, N, D) for N/D close to Q, built as (base rotation) x (a small residual)"""
    k = len(C)
    I = np.eye(k)
    Qb = np.array([[float(int(v)) for v in row] for row in Nb]) / float(Db)
    R = Qb.T @ Q                                   # the residual, close to the identity
    try:
        S = np.linalg.solve((I + R).T, (I - R).T).T
    except np.linalg.LinAlgError:
        return (0, None, None)
    KR = 0.5 * (S - S.T)
    if not np.isfinite(KR).all():
        return (0, None, None)
    thr = np.sqrt(3.0) / 2.0 * m
    # SCORE ONLY THE POLES THAT COULD CHANGE.  The residual is small, so a pole that clears
    # the base rotation by a wide margin clears every candidate here too; ranking is decided
    # by the few near the threshold.  Scoring all of them was 29768 x 27720 per denominator
    # at k = 21 -- ten seconds each, four minutes a pass -- for a stage meant to recover a few
    # dozen poles.  Nothing exact rests on this: the ranking picks which candidates to test,
    # and cleared() tests them over the WHOLE source.
    blk = max(1, int(8e6 // max(1, len(Zf))))
    base = np.empty(len(Tf))
    TB = Tf @ Qb.T
    for tt in range(0, len(Tf), blk):
        base[tt:tt + blk] = np.abs(TB[tt:tt + blk] @ Zf.T).max(1)
    risk = np.nonzero(base > 0.90 * thr)[0]
    safe = len(Tf) - len(risk)
    Tr = Tf[risk]
    nT = len(Tr)
    pool = []
    for d in (denoms or DENOMS):
        if d * np.abs(KR).max() > 9e15:
            continue
        Ki = np.triu(np.rint(d * KR).astype(np.int64), 1)
        N0, D0 = cayley((Ki - Ki.T).tolist(), (d * C).astype(np.int64))
        if N0 is None:
            continue
        N = np.array(Nb, dtype=object) @ np.array(N0, dtype=object)
        D = int(Db) * int(D0)
        # PYTHON integers throughout: N is an object array whose entries run past int64
        # once two rotations are composed, and np.gcd refuses to convert them.
        g = int(D)
        for row in N:
            for v in row:
                g = gcd(g, int(abs(v)))
                if g == 1:
                    break
            if g == 1:
                break
        if g > 1:
            N = N // int(g)
            D //= int(g)
        Qf = np.array([[float(int(v)) for v in row] for row in N]) / float(D)
        TQ = Tr @ Qf.T
        sc = safe + sum(int((np.abs(TQ[t:t + blk] @ Zf.T).max(1) <= thr).sum())
                        for t in range(0, nT, blk))
        pool.append((sc, D, N))
    pool.sort(key=lambda t: (-t[0], t[1]))
    best = (0, None, None)
    for _sc, D, N in pool[:keep]:
        try:
            cnt = exact(N, D)
        except (AssertionError, OverflowError):
            continue
        if cnt > best[0] or (cnt == best[0] and cnt and D < best[2]):
            best = (cnt, N, D)
    return best


def shipped(k):
    """the rotation axis_rotate.py already found for the CAP set, as a starting point"""
    p = os.path.join(DATA, 'poles_rot_%d.npz' % k)
    if not os.path.exists(p):
        return None, None
    z = np.load(p, allow_pickle=True)
    return np.array(z['N'].tolist(), dtype=object), int(z['D'][0])


def improve(k, T, V, starts=4, seed=0, steps=900):
    """(count, N, D, keep) -- the best exact layer this file can build"""
    n = V.shape[1]
    m = int(V[0] @ V[0])
    Tn, Zn = unit(T.astype(np.float64)), unit(V.astype(np.float64))
    # rationalise() ranks candidates in floating point against a single threshold, so hand it
    # the source rescaled to one length; the exact test below keeps the two norms apart.
    Tf = T.astype(np.float64) * np.sqrt(m / (T * T).sum(1))[:, None]
    Vf = V.astype(np.float64)
    C = np.eye(n, dtype=np.int64)

    def exact(N, D):
        return int(cleared(T, V, N, D).sum())

    best = (0, None, None)
    N0, D0 = shipped(k)
    if N0 is not None and N0.shape[0] != n:
        # k = 21's shipped rotation is a 24 x 24 matrix of the AMBIENT coordinates, built by
        # scripts/axis_block.py because no isometry of R^21 was available in that family.
        # Here the frame makes R^21 itself the working space, so that matrix is the wrong
        # shape and there is no rotation to start from.
        N0, D0 = None, None
    if N0 is not None:
        c = exact(N0, D0)
        print('      the shipped rotation (D = %d) carries %d of %d' % (D0, c, len(T)))
        best = (c, N0, D0)
    ident = np.eye(n, dtype=object)
    c = exact(ident, 1)
    print('      no rotation at all carries %d of %d' % (c, len(T)))
    if c > best[0]:
        best = (c, ident, 1)
    rng = np.random.default_rng(seed)
    for t in range(starts):
        if best[0] == len(T):
            break
        if t == 0 and N0 is not None:
            Q = np.array([[float(int(v)) for v in row] for row in N0]) / D0
        else:
            # DETERMINANT +1.  A QR factor is orthogonal but half of them are reflections,
            # and the Cayley transform represents rotations only -- rationalise() would then
            # be approximating something outside its chart, and the start is wasted.
            Q = np.linalg.qr(rng.normal(size=(n, n)))[0]
            Q = Q * np.sign(np.linalg.det(Q))
        t0 = time.time()
        Q = minimax(Tn, Zn, Q, steps=steps)
        w = 0.0
        for s in range(0, len(Tn), 2000):
            w = max(w, float(np.abs((Tn[s:s + 2000] @ Q.T) @ Zn.T).max()))
        cnt, N, D = rationalise(Q, C, Tf, Vf, m, exact, rng=rng, top=3, keep=4)
        # Refine against a rotation that is ALREADY rational, then against the result, and
        # so on: each pass leaves a smaller residual for the next rounding to represent.
        # The denominator multiplies every pass, so this stops on its own -- cleared()
        # rejects a rotation the shipped verifier could not read and refine() drops it.
        for Nb, Db in ([(N, D)] + ([(N0, D0)] if N0 is not None else [])
                       + ([(best[1], best[2])] if best[1] is not None else [])):
            if Nb is None:
                continue
            cb, Nc, Dc = cnt, N, D
            for _pass in range(3):
                c2, N2, D2 = refine(Q, Nb if _pass == 0 else Nc, Db if _pass == 0 else Dc,
                                    C, exact, Tf, Vf, m)
                if c2 <= cb:
                    break
                cb, Nc, Dc = c2, N2, D2
                if cb == len(T):
                    break
            if cb > cnt:
                cnt, N, D = cb, Nc, Dc
            if cnt == len(T):
                break
        print('      start %d: worst cosine %.6f (%.4f of the threshold), rational layer %d'
              '  (%.0f s)' % (t, w, w / THR, cnt, time.time() - t0))
        if cnt > best[0] or (cnt == best[0] and cnt and D < best[2]):
            best = (cnt, N, D)
    keep = np.nonzero(cleared(T, V, best[1], best[2]))[0].astype(np.int32)
    assert len(keep) == best[0]
    return best[0], best[1], best[2], keep


def main():
    argv = sys.argv[1:]
    write = '--write' in argv
    ks = [int(a) for a in argv if a.isdigit()] or list(KS)
    rows = {r['k']: r for r in json.load(open(os.path.join(PKG, 'rows81-95.json')))}
    print(__doc__)
    print('   k  dim   source   shipped   this run   at tau(k)?')
    for k in ks:
        Zk, T, split, code, M, mu, oct_, n = source(k)
        assert is_code(T, split), 'the source is not a 60-degree code'
        cnt, N, D, keep = improve(k, T, Zk)
        old = int(rows[k]['poles'])
        tau = len(T)
        print('  %2d  %3d  %6d  %8d  %9d  %s'
              % (k, 72 + k, len(T), old, cnt, 'YES' if cnt == tau else '%d short' % (tau - cnt)))
        if write and cnt > old:
            np.savez(os.path.join(DATA, 'poles_rec_%d.npz' % k),
                     N=np.array(N.tolist(), dtype=np.int64), D=np.array([D], dtype=np.int64),
                     keep=keep, code=np.array(code, dtype=np.int64),
                     M=np.array(M, dtype=np.int64))
            print('      wrote data/poles_rec_%d.npz  (%d poles, D = %d, frame norm %d)'
                  % (k, cnt, D, mu))


if __name__ == '__main__':
    main()
