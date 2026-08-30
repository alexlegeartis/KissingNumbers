#!/usr/bin/env python3
"""Exhaustive exact verification: every pair, no structural assumptions at all.

    python scripts/verify_exhaustive.py [data/dimension27_200540.txt]

scripts/verify_configuration.py already proves the theorem, quickly and without
third-party packages, but it does so through the geometry: it certifies that the
heads are Leech minimal vectors and then leans on the classical fact that the
Leech lattice has minimal norm 4.  This script assumes none of that.  It knows
only that the coordinates lie in Z[sqrt2,sqrt3,sqrt6], and it compares all
20108045530 pairs.

The two scripts are deliberately different arguments for the same statement.  If
they disagree, something is wrong.

Method.  Write a row as v = (y, tau), where y collects the columns that are
rational in every row and tau the rest, and let n_v be the squared norm.  Then

    <v/sqrt(n_v), w/sqrt(n_w)>  <=  1/2
        <=>  2 (y.y' + tau.tau')  <=  sqrt(n_v n_w)
        <=>  y.y'  <=  floor( sqrt(n_v n_w)/2 - tau.tau' ),

and the right-hand side depends only on the pair of classes (n_v, tau), of which
there are few.  Each threshold is computed exactly by common.floor_half_sqrt_minus.
The remaining comparison is between integers, and the integer dot products are
evaluated with a float64 matrix product, which is exact here because every
partial sum is an integer far below 2^53; the script derives that bound from the
data and checks the product against exact integer arithmetic before using it.

This is the only script in the package that needs numpy.
"""
import collections
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (Report, parse_configuration, sha256sum, rdot, is_rational,
                    floor_half_sqrt_minus)

try:
    import numpy as np
except ImportError:                                          # pragma: no cover
    sys.exit("verify_exhaustive.py needs numpy; "
             "scripts/verify_configuration.py does not, and proves the same thing.")

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, os.pardir, "data", "dimension27_200540.txt")
BA, BB = 4096, 8192


def main(path):
    t0 = time.time()
    rep = Report(f"exhaustive verification of {os.path.basename(path)}")

    rep.section("file")
    rep.note(f"sha256 {sha256sum(path)}")
    cfg = parse_configuration(path)
    n = len(cfg)
    dim = cfg.dimension
    rep.check(cfg.coefficients == [1] * dim,
              "inner product coefficients are all 1 (standard inner product)")
    rep.note(f"dimension {dim}, {n} points, {len(cfg.tokens)} distinct tokens")

    rep.section("columns and classes")
    rational_col = [all(is_rational(row[j]) for row in cfg.rows) for j in range(dim)]
    rat = [j for j in range(dim) if rational_col[j]]
    irr = [j for j in range(dim) if not rational_col[j]]
    rep.check(len(rat) + len(irr) == dim and len(rat) > 0,
              f"{len(rat)} rational columns, {len(irr)} carrying radicals")

    Y = [tuple(row[j][0] for j in rat) for row in cfg.rows]
    TAU = [tuple(row[j] for j in irr) for row in cfg.rows]
    norms = []
    for i in range(n):
        t = rdot(TAU[i], TAU[i])
        if not is_rational(t):
            raise SystemExit(f"row {i}: irrational squared norm")
        norms.append(sum(x * x for x in Y[i]) + t[0])
    rep.check(all(v > 0 for v in norms), "every row is nonzero")
    rep.note(f"{len(set(norms))} distinct squared row norms: {sorted(set(norms))}")

    key = [(norms[i], TAU[i]) for i in range(n)]
    groups = collections.OrderedDict()
    for i, k in enumerate(key):
        groups.setdefault(k, []).append(i)
    gk = list(groups)
    rep.note(f"{len(gk)} classes by (squared norm, radical columns); "
             f"largest has {max(len(v) for v in groups.values())} rows")

    rep.section("exact thresholds")
    G = len(gk)
    THR = [[0] * G for _ in range(G)]
    for a in range(G):
        na, ta = gk[a]
        for b in range(G):
            nb, tb = gk[b]
            THR[a][b] = floor_half_sqrt_minus(na * nb, rdot(ta, tb))
    flat = sorted({THR[a][b] for a in range(G) for b in range(G)})
    rep.note(f"{G}x{G} integer thresholds computed exactly; "
             f"values {flat if len(flat) <= 8 else str(flat[:8]) + ' ...'}")
    rep.note("threshold = floor( sqrt(n_v n_w)/2 - tau.tau' ), an exact algebraic number")

    rep.section("float64 exactness budget")
    mx = max(abs(v) for row in Y for v in row)
    bound = len(rat) * mx * mx
    rep.check(bound < 2 ** 53,
              f"|coordinate| <= {mx}, so every partial sum is an integer <= {bound}",
              f"2^53 = {2 ** 53}")
    order = [i for k in gk for i in groups[k]]
    gof = [0] * n
    for a, k in enumerate(gk):
        for i in groups[k]:
            gof[i] = a
    Ys = np.array([Y[i] for i in order], dtype=np.float64)
    gs = np.array([gof[i] for i in order], dtype=np.int64)
    Yi = np.array([Y[i] for i in order], dtype=np.int64)
    rng = np.random.default_rng(0)
    ia, ib = rng.integers(0, n, 512), rng.integers(0, n, 512)
    rep.check((Ys[ia] @ Ys[ib].T == (Yi[ia] @ Yi[ib].T).astype(np.float64)).all(),
              "float64 matrix product is bit-exact against int64 on a 512x512 sample")

    rep.section("the sweep")
    THRa = np.array(THR, dtype=np.float64)
    bad, worst, t1 = 0, -10 ** 18, time.time()
    for a in range(0, n, BA):
        Ab, ga = Ys[a:a + BA], gs[a:a + BA]
        for b in range(a, n, BB):
            Bb, gb = Ys[b:b + BB], gs[b:b + BB]
            g = Ab @ Bb.T
            # Blank the diagonal wherever the two blocks overlap: row a+i of the
            # A block is column a+i-b of the B block.  The column mask decides
            # this on its own, so it stays correct for any BA and BB -- an outer
            # guard such as "b <= a < b + len(Bb)" would only be right while
            # BB >= BA keeps the blocks from overlapping partially.
            i = np.arange(len(Ab))
            c = a + i - b
            ok = (c >= 0) & (c < g.shape[1])
            g[i[ok], c[ok]] = -10 ** 18
            # subtract in place: g is up to BA x BB float64, and a separate
            # difference array would be a third buffer of that size
            g -= THRa[np.ix_(ga, gb)]
            bad += int((g > 0).sum())
            worst = max(worst, float(g.max()))
        if a % (BA * 12) == 0:
            print(f"         {a}/{n}  violations {bad}  worst margin {int(worst)}  "
                  f"{time.time() - t1:.0f}s", flush=True)
    rep.check(bad == 0, f"all {n * (n - 1) // 2} pairs satisfy <p,q> <= 1/2",
              f"max (y.y' - threshold) = {int(worst)}, so the bound is attained"
              if worst == 0 else f"max margin {int(worst)}")

    print()
    print(f"  elapsed {time.time() - t0:.0f}s")
    return rep.done(f"VERIFIED EXHAUSTIVELY: {n} distinct unit vectors in R^{dim},\n"
                    f"           pairwise inner product at most 1/2, so K({dim}) >= {n}.")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT))
