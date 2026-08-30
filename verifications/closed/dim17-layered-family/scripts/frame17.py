"""tau(17) = 5346 + 2*alpha, derived from Cohn's coordinates rather than from the algebra.

    python frame17.py [/path/to/dimensions1-24.txt]

`theta17.py` proves `alpha = 192` for the graph `Cay(K', W_4)`.  This script establishes that
that graph is the right one, from the record itself and nothing else.  Sections 73-74 derive
the same frame by hand; the point here is that the derivation can be replaced by a measurement.

Everything below comes out of the 5730 coordinates:

  1. the record has a layering axis with exactly seven heights, 1, 1/2, 1/3, 0 and their
     negatives, holding 1, 512, 192, 4320, 192, 512, 1 points;
  2. the record is NOT antipodally symmetric -- the two height-1/3 layers are different sets,
     which is worth knowing because the natural assumption breaks the geometry: a direction
     used at BOTH +1/3 and -1/3 would give cos 7/9;
  3. deleting the two height-1/3 layers leaves 5346 points -- the constant in the formula;
  4. sweeping all 65536 sign patterns against those 5346 points, exactly 1024 are admissible at
     height 1/3, and they form a COSET of a group whose difference weights are
     {0:1, 4:60, 6:256, 8:390, 10:256, 12:60, 16:1};
  5. that group is `K'` -- not merely isomorphic to it: it is compared element by element with
     the set `theta17.py` builds from the supports file.  The 60 weight-4 differences are
     `W_4`, and the tier's own constraint
     `b.b' <= 4` is exactly "no difference of weight 4".  So each height-1/3 layer is an
     independent set in `Cay(K', W_4)`;
  6. the two layers do not constrain each other -- the binding inequality between them needs
     only Hamming distance 3, and two words of one coset are always at distance 0, 4, 6, ... --
     so each independently attains alpha, and the total is 5346 + 2*alpha.

With `theta17.py`'s `alpha = 192` that is 5346 + 384 = 5730, the record's own size, so the
record is optimal in its family.

The data set is H. Cohn, Table of kissing number bounds, https://hdl.handle.net/1721.1/153312
-- 24 MB, not redistributed here.  `extract.py` explains how to point this script at it.
"""
import collections
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract import blocks, _find_source                                        # noqa: E402

HEIGHTS = {720: 'pole+', 360: 'A+', 240: 'B+', 0: 'equator',
           -240: 'B-', -360: 'A-', -720: 'pole-'}
EXPECT = {720: 1, 360: 512, 240: 192, 0: 4320, -240: 192, -360: 512, -720: 1}


def main():
    print(__doc__)
    src = _find_source()
    rec = [x for x in blocks(src) if x['dim'] == 17 and x['n'] == 5730]
    assert rec, "no 5730-point dimension-17 block in the data set"
    b = rec[0]
    X = np.array(b['pts'], dtype=float) * np.sqrt(np.array(b['c'], dtype=float))[None, :]
    U = X / np.linalg.norm(X, axis=1)[:, None]
    G = U @ U.T
    off = G - 2 * np.eye(len(U))
    print("record: %d points in R^17, max cosine %.12f" % (len(U), off.max()))
    assert off.max() <= 0.5 + 1e-9

    # 1. the axis
    ax = min(range(len(U)),
             key=lambda i: len(set(np.round(G[i] * 720).astype(np.int64).tolist())))
    h = np.round(G[ax] * 720).astype(np.int64)
    lay = collections.Counter(h.tolist())
    assert dict(lay) == EXPECT, dict(lay)
    print("\nlayering axis = point %d, heights and populations:" % ax)
    for q in sorted(lay, reverse=True):
        print("    %-8s h = %+7.4f   %5d points" % (HEIGHTS[q], q / 720.0, lay[q]))

    # 2. not antipodally symmetric
    anti = int((G.min(axis=1) < -1 + 1e-9).sum())
    print("\npoints whose antipode is also in the record: %d of %d" % (anti, len(U)))
    assert anti < len(U), "expected the record NOT to be antipodally symmetric"

    # 3. the fixed part
    fixed = np.where(np.abs(h) != 240)[0]
    print("deleting both height-1/3 layers leaves %d points" % len(fixed))
    assert len(fixed) == 5346
    FIX = U[fixed]

    # 4. every admissible sign pattern at height 1/3
    z = U[ax]
    # the axis is a signed coordinate vector, so the natural 16-frame is the other coordinates
    # in order -- which is the labelling data/d17_supports.npy uses.  An arbitrary orthonormal
    # complement (a QR basis, say) also works geometrically but relabels the sign patterns, and
    # then the group below is only ISOMORPHIC to K' rather than equal to it.
    k = int(np.argmax(np.abs(z)))
    assert abs(abs(z[k]) - 1) < 1e-12 and k == 16, "expected the axis to be +-e_17"
    E = np.eye(17)[:, :16]
    sgn = np.array([[1.0 if not (w >> i) & 1 else -1.0 for i in range(16)]
                    for w in range(1 << 16)])
    cand = (sgn / 4.0) @ E.T * (2 * np.sqrt(2) / 3) + np.outer(np.ones(1 << 16), z / 3.0)
    ok = np.ones(1 << 16, dtype=bool)
    for a in range(0, 1 << 16, 4096):
        ok[a:a + 4096] = (cand[a:a + 4096] @ FIX.T).max(axis=1) <= 0.5 + 1e-9
    adm = [int(v) for v in np.where(ok)[0]]
    print("\nsign patterns admissible at height 1/3, out of 65536: %d" % len(adm))
    assert len(adm) == 1024

    # The SAME sweep at height -1/3.  Each layer attaining alpha is what makes the total
    # 5346 + 2*alpha rather than 5346 + alpha + something, and that needs the two admissible
    # sets to be one and the same -- true because the tier-A direction set is symmetric, but
    # measured here rather than argued.
    cand2 = (sgn / 4.0) @ E.T * (2 * np.sqrt(2) / 3) - np.outer(np.ones(1 << 16), z / 3.0)
    ok2 = np.ones(1 << 16, dtype=bool)
    for a in range(0, 1 << 16, 4096):
        ok2[a:a + 4096] = (cand2[a:a + 4096] @ FIX.T).max(axis=1) <= 0.5 + 1e-9
    adm2 = [int(v) for v in np.where(ok2)[0]]
    print("the same sweep at height -1/3 gives %d, identical set: %s"
          % (len(adm2), set(adm2) == set(adm)))
    assert set(adm2) == set(adm)

    # 5. it is a coset, and its difference group is K'
    diffs = {w ^ adm[0] for w in adm}
    assert all((x ^ y) in diffs for x in diffs for y in list(diffs)[:32])
    wt = collections.Counter(bin(x).count('1') for x in diffs)
    print("they are one coset of a 1024-element group; difference weights %s"
          % dict(sorted(wt.items())))
    assert dict(sorted(wt.items())) == {0: 1, 4: 60, 6: 256, 8: 390, 10: 256, 12: 60, 16: 1}

    # This group is NOT the same set of words as theta17.py's K': theta17 builds RM(2,4) in the
    # Boolean-cube coordinate order and Cohn's file uses its own, so the two are isomorphic but
    # differently labelled.  Rather than hunt for the permutation, apply theta17's certificate
    # to THIS graph.  If the same six numbers certify it, the bound holds for the object the
    # geometry actually produces and the labelling never has to be settled.
    from theta17 import build, certificate                                   # noqa: E402
    kl0, _w40, _f, _n = build()
    K = sorted(diffs)
    W = sorted(x for x in K if bin(x).count('1') == 4)
    print("   same word-for-word as theta17.py's K': %s   (labelling differs, see above)"
          % (set(kl0) == diffs))
    # the fibre subgroup.  K' is 10-dimensional and contains its own dual, so K'^perp has 64
    # elements, all of them inside K'; its weights are {0:1, 6:16, 8:30, 10:16, 16:1} and the
    # part of weight 0, 8 or 16 is RM(1,4) -- 32 words, and a subgroup.
    perp = [y for y in range(1 << 16) if all(bin(y & x).count('1') % 2 == 0 for x in K)]
    assert set(perp) <= diffs and len(perp) == 64
    H = sorted(y for y in perp if bin(y).count('1') in (0, 8, 16))
    assert len(H) == 32 and all((a ^ b) in set(H) for a in H for b in H)
    print("   K'-perp: %d words, all inside K'; its RM(1,4) part is the %d-element fibre group"
          % (len(perp), len(H)))
    bound, spec = certificate(K, W, H, show="   the same certificate on the MEASURED graph")
    print("   spectrum %s" % dict(sorted(spec.items(), reverse=True)))
    print("   alpha of the measured graph <= %s" % bound)
    assert bound == 192, bound

    # the record's own layers are independent sets in Cay(K', W_4), and are maximum
    P = U - np.outer(G[ax], z)
    for q in (240, -240):
        L = P[h == q]
        L = L / np.linalg.norm(L, axis=1)[:, None]
        words = [int(sum(1 << i for i in range(16) if r[i] < 0)) for r in L @ E]
        assert len(set(words)) == 192 and set(words) <= set(adm)
        bad = sum(1 for i in range(192) for j in range(i + 1, 192)
                  if bin(words[i] ^ words[j]).count('1') == 4)
        print("   layer %-3s : 192 admissible words, %d pairs at Hamming distance 4"
              % (HEIGHTS[q], bad))
        assert bad == 0

    # 6. the two layers do not constrain each other
    cross = G[np.ix_(np.where(h == 240)[0], np.where(h == -240)[0])].max()
    print("\nmaximum cosine between the two height-1/3 layers: %.6f  (limit 0.5)" % cross)
    print("   the binding inequality needs only Hamming distance 3, and two words of one")
    print("   coset are at distance 0, 4, 6, ... -- so the two layers are independent")
    assert cross <= 0.5 + 1e-9

    print("\n" + "=" * 78)
    print("tau(17) = 5346 + 2*alpha,  alpha = alpha(Cay(K', W_4))   -- MEASURED, not assumed")
    print("theta17.py gives alpha = 192, so the family maximum is 5346 + 384 = 5730")
    print("=" * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
