#!/usr/bin/env python3
"""Parse the Catalogue of Lattices pages for the extremal even unimodular 48-dimensional
lattices P_48p and P_48q, and verify them in exact integer arithmetic.

    http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/P48p.html   (SL2(23)xS3):2
    http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/P48q.html   SL2(47)

(There is no P48n.html on the server; a 404 is returned.  P_48p and P_48q are the two
classical extremal lattices of SPLAG p. 195.)

Checks performed, all in Python integers (no floating point):
  * Gram is 48x48, symmetric, even diagonal (so the lattice is even);
  * determinant EXACTLY 1, by fraction-free Bareiss elimination;
  * every published group generator A satisfies A G A^T = G exactly, det A = +-1.

Writes  p48<x>_gram.npy  (48x48 int64)  and  p48<x>_gens.npy  (g x 48 x 48 int64).

    python parse_p48.py <page.html> <tag> <outdir>
"""
import re, sys, os
import numpy as np


def strip_html(s):
    return re.sub(r'<[^>]*>', ' ', s)


def ints_between(raw, a, b):
    blk = strip_html(raw[raw.index(a):raw.index(b)])
    return [int(t) for t in blk.split() if re.fullmatch(r'-?\d+', t)]


def bareiss_det(M):
    """exact determinant of an integer matrix by fraction-free Bareiss elimination"""
    A = [[int(x) for x in row] for row in M]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = A[i][j] * A[k][k] - A[i][k] * A[k][j]
                assert num % prev == 0, "Bareiss division not exact"
                A[i][j] = num // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def main():
    if len(sys.argv) < 4:
        raise SystemExit(
            "usage: python parse_p48.py <page.html> <tag> <outdir>\n"
            "\n"
            "  e.g. python parse_p48.py ../data/P48p.html p48p ../data\n"
            "\n"
            "Parses a Catalogue of Lattices page into a Gram matrix and generator matrices\n"
            "and verifies them exactly: symmetric, even diagonal, determinant 1 by\n"
            "fraction-free Bareiss elimination, and A G A^T = G for every generator.\n"
            "\n"
            "The outputs are already shipped in ../data/; this script is here so that the\n"
            "provenance can be re-checked from the catalogue page itself.")
    src, tag, out = sys.argv[1], sys.argv[2], sys.argv[3]
    raw = open(src, 'r', encoding='latin-1').read()

    nums = ints_between(raw, '<a NAME="GRAM">', '<a NAME="MINIMAL_NORM">')
    # two catalogue layouts: "48 48" + full square (P_48p), or "48" + lower triangle (P_48q)
    if nums[0] == 48 and nums[1] == 48 and len(nums) == 2 + 48 * 48:
        G = np.array(nums[2:], dtype=np.int64).reshape(48, 48)
    else:
        assert nums[0] == 48, nums[:4]
        tri = nums[1:]
        assert len(tri) == 48 * 49 // 2, len(tri)
        G = np.zeros((48, 48), dtype=np.int64)
        k = 0
        for i in range(48):
            for j in range(i + 1):
                G[i, j] = G[j, i] = tri[k]
                k += 1

    nums = ints_between(raw, '<a NAME="GROUP_GENERATORS">', '<a NAME="PROPERTIES">')
    ngen = nums[0]
    rest = nums[1:]
    assert len(rest) == ngen * (2 + 48 * 48), (len(rest), ngen)
    GENS = np.empty((ngen, 48, 48), dtype=np.int64)
    for g in range(ngen):
        base = g * (2 + 48 * 48)
        assert rest[base] == 48 and rest[base + 1] == 48
        GENS[g] = np.array(rest[base + 2:base + 2 + 48 * 48],
                           dtype=np.int64).reshape(48, 48)

    # ---------- exact verification ----------
    assert G.shape == (48, 48)
    assert (G == G.T).all(), "Gram not symmetric"
    d = np.diag(G)
    assert (d % 2 == 0).all(), "Gram diagonal not even"
    det = bareiss_det(G)
    print("%s: 48x48 symmetric, diagonal %s, Bareiss determinant = %d"
          % (tag, sorted(set(int(x) for x in d)), det))
    assert det == 1, "determinant is not 1"
    for i, A in enumerate(GENS):
        assert (A @ G @ A.T == G).all(), "generator %d is not an isometry" % i
        da = bareiss_det(A)
        assert abs(da) == 1
        print("   generator %d: A G A^T == G exactly, det A = %d, max|A_ij| = %d"
              % (i, da, int(np.abs(A).max())))

    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, '%s_gram.npy' % tag), G)
    np.save(os.path.join(out, '%s_gens.npy' % tag), GENS)
    print("   saved %s_gram.npy %s_gens.npy  (%d generators)" % (tag, tag, ngen))


if __name__ == '__main__':
    main()
