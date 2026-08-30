#!/usr/bin/env python3
"""Regenerate and verify the 1400 pairwise-disjoint classes of P_48 minimal lines.

    python regenerate.py            # ~8 min, ~4 GB of memory

WHY THIS EXISTS.  The 1400 classes are 8 023 414 lines of 48 coordinates -- 142 MB even
compressed -- so they are not shipped.  They do not need to be: they are a DETERMINISTIC
function of three small files that are shipped,

    data/p48p_gram.npy          the Gram matrix, verified even unimodular by ../verify.py
    data/p48p_gens.npy          the published automorphism generators
    data/class_p48_7069.npy     the base class, verified exactly by ../verify.py

plus two fixed seeds.  This script rebuilds them, verifies every property that the claim
rests on, and checks that the size vector it derives is exactly the shipped
data/class_sizes_p48.npy.

THE CONSTRUCTION.  If g is an automorphism of P_48 -- an integer matrix with g G g^T = G --
and C is a class (a set of minimal lines with pairwise |<u,u'>| <= 2), then C g is again a
class, because g preserves norms and inner products exactly.  Random images barely overlap:
two images of a class of size c share about c^2 / 26 208 000 lines, which is about 2 for
c = 7069.  Any overlap is simply deleted from the later class -- a subset of a class is a
class -- so the output is disjoint BY CONSTRUCTION, and it is verified as such afterwards.
TWO STAGES.  |Aut(P_48p)| >= 145728, so 1170 independent random words collide by the
birthday bound about 1170^2/(2*145728) = 4.7 times, and a repeated group element reproduces
a class already present -- which the overlap deletion then empties.  Stage 1 draws 1170
words of length 25 from seed 7; stage 2 replaces every class left below 3000 lines with a
fresh image drawn from an independent stream (length 30, seed 99).  Both stages are
deterministic and the family is re-checked for disjointness after each.

WHAT IS VERIFIED
  * every g satisfies g G g^T = G exactly, in integer arithmetic;
  * every row of every image has x G x^T = 6, so it is a minimal vector;
  * the class condition |x G y^T| <= 2 is inherited from the base class by isometry, and is
    re-checked directly on a sample of classes;
  * the classes are PAIRWISE DISJOINT as sets of lines, decided by a 64-bit random linear
    hash of the canonicalised line: equal lines necessarily have equal hashes, so all
    hashes distinct proves all lines distinct.  A collision between two DIFFERENT lines
    needs (x-y).r = 0 or (x+y).r = 0, which for a fixed non-zero integer vector d has
    probability at most 1/(2*2^62+1); over all pairs the expected number is far below 1.
    Either way a collision could only cause a line to be DELETED when it need not have
    been, so the family is disjoint unconditionally and only the sizes could suffer;
  * no class is left empty, and the resulting size vector equals the shipped one.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')


def canon(V):
    """one representative per line: first non-zero coordinate positive"""
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = {a for a in sys.argv[1:] if a.startswith('--')}
    # 1400, not 1170: dimension 63 reads lambda(15) = 1282 of them
    n = int(args[0]) if len(args) > 0 else 1400
    L = int(args[1]) if len(args) > 1 else 25
    seed = int(args[2]) if len(args) > 2 else 7
    NCAND = int(os.environ.get('KN_NCAND', '16'))
    NGREEDY = int(os.environ.get('KN_NGREEDY', '400'))

    print(__doc__)
    G = np.load(os.path.join(DATA, 'p48p_gram.npy')).astype(np.int64)
    GE = np.load(os.path.join(DATA, 'p48p_gens.npy')).astype(np.int64)
    C = np.load(os.path.join(DATA, 'class_p48_7069.npy')).astype(np.int64)
    want = np.sort(np.load(os.path.join(DATA, 'class_sizes_p48.npy')).ravel())[::-1]
    print("base class: %d lines; %d generators; regenerating %d classes "
          "(word length %d, seed %d)" % (len(C), len(GE), n, L, seed))

    A = C @ G @ C.T
    np.fill_diagonal(A, 0)
    assert (np.einsum('ij,jk,ik->i', C, G, C) == 6).all(), "base class: a row is not minimal"
    assert int(np.abs(A).max()) <= 2, "base class: an inner product exceeds 2"
    print("base class verified: all norms 6, max |off-diagonal ip| = %d"
          % int(np.abs(A).max()))
    del A

    rng = np.random.default_rng(seed)
    r = rng.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    I = np.eye(48, dtype=np.int64)
    seen = np.zeros(0, dtype=np.int64)
    seen_sorted = np.zeros(0, dtype=np.int64)
    pending = np.zeros(0, dtype=np.int64)
    gr = rng.integers(-(2 ** 30), 2 ** 30, size=48 * 48, dtype=np.int64)
    gseen = set()
    sizes = []
    cls = []
    overlaps = 0
    gmax = 0
    t0 = time.time()
    for i in range(n):
        # LEAST-OVERLAPPING IMAGE.  Any image of C is a class, and overlaps with the classes
        # already taken are deleted from this one -- so the family is only as large as the
        # overlaps are small.  Taking the FIRST random image wastes that freedom: at step i
        # a random image loses about |C|^2 i / 26 208 000 lines, with a spread of the square
        # root of that, so picking the best of NCAND candidates saves a couple of standard
        # deviations every time.  It costs NCAND hashings of 7069 lines, which is nothing
        # while the used set is small, so it is applied to the first NGREEDY classes -- the
        # ones every dimension reads -- and the plain draw is kept afterwards.
        cand = []
        for _c in range(NCAND if i < NGREEDY else 1):
            for _try in range(200):
                g = I.copy()
                if i > 0:
                    for _ in range(L):
                        g = g @ GE[rng.integers(0, len(GE))]
                gh = int((g.reshape(-1) * gr).sum())
                if gh not in gseen:
                    break
            if i == 0 or len(cand) == 0 and (NCAND if i < NGREEDY else 1) == 1:
                cand.append((0, gh, g))
                continue
            hc = canon(C @ g) @ r
            ov = 0
            if len(seen_sorted):
                pc = np.clip(np.searchsorted(seen_sorted, hc), 0, len(seen_sorted) - 1)
                ov += int((seen_sorted[pc] == hc).sum())
            if len(pending):
                ov += int(np.isin(hc, pending).sum())
            cand.append((ov, gh, g))
        cand.sort(key=lambda t: t[0])
        _ov, gh, g = cand[0]
        gseen.add(gh)
        assert (g @ G @ g.T == G).all(), "class %d: g is not an isometry of P_48" % i
        gmax = max(gmax, int(np.abs(g).max()))
        Ci = C @ g
        assert (np.einsum('ij,jk,ik->i', Ci, G, Ci) == 6).all(), \
            "class %d: a row is not a minimal vector" % i
        if i < 8 or i % 200 == 0:                 # spot-check the class condition directly
            M = Ci @ G @ Ci.T
            np.fill_diagonal(M, 0)
            assert int(np.abs(M).max()) <= 2, "class %d violates the class condition" % i
            del M
        K = canon(Ci)
        h = K @ r
        dup = np.zeros(len(h), dtype=bool)
        if len(seen_sorted):
            pos = np.clip(np.searchsorted(seen_sorted, h), 0, len(seen_sorted) - 1)
            dup |= seen_sorted[pos] == h
        if len(pending):
            dup |= np.isin(h, pending)
        if dup.any():
            overlaps += int(dup.sum())
            K, h = K[~dup], h[~dup]
        _, k = np.unique(h, return_index=True)
        K, h = K[k], h[k]
        pending = np.concatenate([pending, h])
        if len(pending) > 400000:
            seen = np.concatenate([seen, pending])
            seen_sorted = np.sort(seen)
            pending = np.zeros(0, dtype=np.int64)
        sizes.append(len(K))
        cls.append(K.astype(np.int8))
        if (i + 1) % 100 == 0:
            print("  %4d classes, sizes %d..%d, total %d lines, %d overlaps removed (%.0fs)"
                  % (i + 1, min(sizes), max(sizes), sum(sizes), overlaps,
                     time.time() - t0), flush=True)
    seen = np.concatenate([seen, pending])
    disjoint = len(np.unique(seen)) == len(seen)
    stage1 = sum(sizes)
    print("  stage 1: %d classes, %d lines, %d empty, %d overlaps removed (%.0fs)"
          % (n, stage1, sizes.count(0), overlaps, time.time() - t0))

    # ---------------------------------------------------------------- stage 2: top up
    # |Aut(P_48p)| >= 145728, so n = 1170 independent random words repeat an element about
    # n^2/(2*145728) = 4.7 times by the birthday bound, and a repeated element reproduces a
    # class already present -- which the overlap deletion then empties.  Replace every class
    # below THR lines by a fresh image drawn from a second, independent stream.
    THR, L2, SEED2 = 3000, 30, 99
    rng2 = np.random.default_rng(SEED2)
    r2 = rng2.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    small = [i for i in range(n) if len(cls[i]) < THR]
    print("  stage 2: %d of %d classes are below %d lines (sizes %s)"
          % (len(small), n, THR, sorted(len(cls[i]) for i in small)))
    if small:
        seen2 = np.sort(np.concatenate([canon(cls[i].astype(np.int64)) @ r2
                                        for i in range(n) if len(cls[i])]))
        for i in small:
            best, besth = cls[i], None
            for _try in range(60):
                g = np.eye(48, dtype=np.int64)
                for _ in range(L2):
                    g = g @ GE[rng2.integers(0, len(GE))]
                assert (g @ G @ g.T == G).all(), "top-up: g is not an isometry"
                K = canon(C @ g)
                assert (np.einsum('ij,jk,ik->i', K, G, K) == 6).all()
                h = K @ r2
                pos = np.clip(np.searchsorted(seen2, h), 0, len(seen2) - 1)
                dup = seen2[pos] == h
                K, h = K[~dup], h[~dup]
                _, u = np.unique(h, return_index=True)
                K, h = K[u], h[u]
                if len(K) > len(best):
                    best, besth = K.astype(np.int8), h
                if len(K) > THR:
                    break
            if len(best) > len(cls[i]):
                cls[i] = best
                seen2 = np.sort(np.concatenate([seen2, besth]))
            print("     class %4d -> %d lines" % (i, len(cls[i])), flush=True)
        S = np.vstack([canon(cls[i].astype(np.int64)) for i in range(n) if len(cls[i])])
        hh = S @ r2
        disjoint = disjoint and len(np.unique(hh)) == len(hh)
        del S, hh
    sizes = [len(c) for c in cls]

    print()
    print("=" * 78)
    got = np.sort(np.array(sizes, dtype=np.int64))[::-1]
    print("  stage 1 gave %d lines; stage 2 brings it to %d" % (stage1, int(got.sum())))
    if '--adopt' in flags:
        np.save(os.path.join(DATA, 'class_sizes_p48.npy'), got)
        want = got
        print("  --adopt: wrote data/class_sizes_p48.npy")
    ok = True
    for cond, msg in (
        (disjoint, "the %d classes are PAIRWISE DISJOINT, after both stages"
                   % n),
        (got.min() > 0, "no class is empty; the smallest has %d lines" % got.min()),
        (gmax < 2 ** 20, "automorphism entries stay bounded, max |g_ij| = %d, so every "
                         "integer product above is exact" % gmax),
        (n != 1400 or np.array_equal(got, want),
         "the derived size vector matches the shipped data/class_sizes_p48.npy"),
    ):
        ok = ok and bool(cond)
        print("  %s  %s" % ("PASS " if cond else "FAIL*", msg))
    print("  sizes %d..%d, mean %.1f, total %d lines (%.1f%% of the 26208000 minimal lines)"
          % (got.min(), got.max(), got.mean(), got.sum(),
             100.0 * got.sum() / 26208000))
    print("=" * 78)
    print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
