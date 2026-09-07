#!/usr/bin/env python3
"""Exact verification of a kissing configuration in R^27, from the data file alone.

    python scripts/verify_configuration.py [data/dimension27_200540.txt]

This script proves the theorem directly from the coordinate file.  It does not
use the construction in the paper, and it uses no floating-point arithmetic and
no third-party packages.  If it passes, the file is a valid kissing
configuration and the lower bound follows, even if the paper were misstated.

How the check is made exact and fast
------------------------------------
Every coordinate lies in Z[sqrt2, sqrt3, sqrt6].  Splitting a row as v = (y, t)
along the 24 rational columns and the 3 columns that carry radicals, and
rescaling each row to a canonical length, the rows fall into three families:

    equator   (y, 0)   y.y =  32
    cap       (h, t)   h.h = 288, t.t = 144      (h = 3y with y.y = 32)
    axis      (0, a)   a.a =  48

The 196560 vectors y arising as equator heads and cap heads are then certified
to be the minimal vectors of the Leech lattice, from the data:

  * their shape census is 1104 + 97152 + 98304, the census of the minimal
    vectors of Lambda_24;
  * the supports of the 97152 vectors of shape (+-2^8, 0^16) are 759 octads
    which span a [24, 12] binary code whose weight distribution is that of the
    extended binary Golay code, and which are exactly its weight-8 words;
  * every one of the 196560 satisfies the congruences defining Lambda_24 with
    respect to that code.

Two classical facts are taken as given, and nothing else:

  (a) a binary [24,12] code with weight enumerator 1 + 759x^8 + 2576x^12 +
      759x^16 + x^24 is permutation-equivalent to the extended binary Golay
      code, so the congruences above define a permuted copy of Lambda_24;
  (b) Lambda_24 has minimal norm 4, that is, its minimal vectors have squared
      norm 32 in this scaling (Conway and Sloane, SPLAG chapter 10).

Given those, the set of heads lies in a lattice of minimal norm 32, so for
distinct heads u, v the difference u - v is a nonzero lattice vector,
|u - v|^2 >= 32, and u.v <= 16.  The equator-equator and equator-cap conditions
then need no computation at all.  Everything else reduces to a few hundred
thousand integer dot products inside the four cap classes, and to two 12 x 12
tables in the ring.

Run scripts/verify_exhaustive.py for a check that assumes nothing whatsoever,
including nothing about the Leech lattice; it compares all 20108045530 pairs
directly and needs numpy.
"""
import math
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (Report, parse_configuration, sha256sum, rdot, radd,
                    rstr, shape_str, ring_cmp_rational, ring_cmp,
                    split_families, certify_leech_minimal)

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, os.pardir, "data", "dimension27_200540.txt")


def main(path):
    t0 = time.time()
    rep = Report(f"exact verification of {os.path.basename(path)}")

    # ----------------------------------------------------------------- file
    rep.section("file")
    rep.note(f"sha256 {sha256sum(path)}")
    cfg = parse_configuration(path, expect_dimension=27)
    n = len(cfg)
    rep.check(cfg.dimension == 27, "dimension is 27")
    rep.check(cfg.coefficients == [1] * 27,
              "inner product coefficients are all 1 (standard inner product)")
    rep.note(f"parsed {n} rows, {len(cfg.tokens)} distinct coordinate tokens")

    # ------------------------------------------------- rational columns
    rep.section("coordinate ring and families")
    try:
        equator, caps, axis, rat, irr = split_families(cfg)
    except ValueError as e:
        rep.check(False, "rows split into equator, cap and axis families", str(e))
        return rep.done("")
    rep.check(rat == list(range(24)), "columns 0-23 are rational in every row")
    rep.check(irr == [24, 25, 26], "columns 24-26 are the ones carrying radicals")
    rep.check(len(equator) + len(caps) + len(axis) == n,
              "every row is an equator, cap or axis row, at a canonical scale",
              f"{len(equator)} + {len(caps)} + {len(axis)}")
    rep.check(all(sum(x * x for x in y) == 32 for y in equator),
              "every equator head has squared norm 32")
    rep.check(all(sum(x * x for x in y) == 32 for y, _ in caps),
              "every cap head is 3y with y of squared norm 32")
    rep.check(all(rdot(t, t) == (144, 0, 0, 0) for _, t in caps),
              "every cap tail has squared norm 144")
    rep.check(all(rdot(a, a) == (48, 0, 0, 0) for a in axis),
              "every axis tail has squared norm 48")
    rep.note("so the three families are exactly the unit vectors")
    rep.note("  (y,0)/sqrt(32),  (3y,t)/sqrt(432),  (0,a)/sqrt(48)")

    # ------------------------------------------------------- distinctness
    rep.section("distinctness")
    keys = ([("e",) + y for y in equator]
            + [("c",) + y + tuple(c for tv in t for c in tv) for y, t in caps]
            + [("a",) + tuple(c for tv in a for c in tv) for a in axis])
    rep.check(len(set(keys)) == n,
              f"all {n} points are distinct after rescaling to canonical length",
              f"{len(set(keys))} distinct keys")

    # ------------------------------------------------- the Leech lattice
    rep.section("the 196560 heads are the minimal vectors of the Leech lattice")
    # each class is carried by three directions, so cap heads repeat; dedupe
    distinct_cap_heads = sorted({y for y, _ in caps})
    heads = equator + distinct_cap_heads
    rep.check(len(set(heads)) == len(heads), "the heads are pairwise distinct")
    rep.check(len(heads) == 196560, "there are 196560 distinct heads in all",
              f"{len(equator)} equator + {len(distinct_cap_heads)} distinct cap "
              f"(from {len(caps)} caps)")

    cert = certify_leech_minimal(heads)
    rep.check(cert["norms_32"], "every head has squared norm 32")
    rep.check(cert["census_ok"], "shape census is 1104 + 97152 + 98304",
              "  ".join(f"{v} x {shape_str(k)}" for k, v in sorted(cert["census"].items())))
    rep.check(cert["n_octads"] == 759,
              "the 97152 vectors of shape (+-2^8) have 759 supports", str(cert["n_octads"]))
    rep.check(cert["golay_dim"] == 12,
              "those 759 octads span a 12-dimensional binary code",
              f"|C| = {cert['golay_size']}")
    rep.check(cert["weights_ok"],
              "its weight distribution is that of the extended binary Golay code",
              " + ".join(f"{v}x^{k}" for k, v in cert["weights"].items()))
    rep.check(cert["octads_are_weight8"],
              "the 759 octads are exactly the weight-8 codewords")
    rep.check(cert["all_in_leech"],
              "every head satisfies the congruences defining Lambda_24",
              "" if cert["all_in_leech"] else f"first failure {cert['first_bad']}")
    rep.note("A [24,12] code with that weight enumerator is the extended Golay code up")
    rep.note("to permutation, so these congruences define a permuted Lambda_24, whose")
    rep.note("minimal norm is 4 -- squared norm 32 here.  For distinct heads u, v the")
    rep.note("difference u - v is a nonzero lattice vector, so |u-v|^2 >= 32 and u.v <= 16.")
    LEECH_MAX = 16

    # ---------------------------------------------------------- condition E
    rep.section("condition (E): cap directions are deleted from the equator")
    eq_set = set(equator)
    cap_heads = {y for y, _ in caps}
    rep.check(not (eq_set & cap_heads),
              "no cap head occurs as an equator point",
              f"{len(eq_set)} equator, {len(cap_heads)} distinct cap heads")
    rep.check(len(eq_set) + len(cap_heads) == 196560,
              "equator and cap heads together are all 196560 minimal vectors")

    # ------------------------------------------------ classes and directions
    rep.section("cap directions and their classes")
    by_dir = {}
    for y, t in caps:
        by_dir.setdefault(t, []).append(y)
    dirs = sorted(by_dir)
    rep.check(len(dirs) == 12, "the caps lie on 12 distinct directions")
    sizes = sorted({len(v) for v in by_dir.values()})
    rep.check(sizes == [496], "each direction carries 496 caps")
    classes = []
    class_of = {}
    for d in dirs:
        key = frozenset(by_dir[d])
        if key not in classes:
            classes.append(key)
        class_of[d] = classes.index(key)
    # the number of distinct classes is not prescribed: the case analysis below
    # is general, so this script verifies any configuration of this shape,
    # including the published 200044-point one, which uses five classes.
    rep.check(1 <= len(classes) <= len(dirs),
              f"{len(classes)} distinct classes occur among the {len(dirs)} directions",
              f"assignment {[class_of[d] for d in dirs]}")

    clist = [sorted(c) for c in classes]
    NEG = -10 ** 9
    # Inside a class every pair has to be examined: the Leech certification
    # bounds the dot product by 16, but a shared direction needs 8, so the
    # maximum is genuinely computed here.
    self_max = []
    for ci in clist:
        m = NEG
        for a in range(len(ci)):
            u = ci[a]
            for b in range(a + 1, len(ci)):
                s = sum(map(int.__mul__, u, ci[b]))
                if s > m:
                    m = s
        self_max.append(m)
    # Between two classes the certification above already proves the maximum
    # cannot exceed LEECH_MAX, so the search only has to find a pair attaining
    # it and may stop there.  Either way `m` ends as the true maximum: the loop
    # stops early only at a value that is proved to be the ceiling.
    cross_max, cross_attained = {}, True
    for i in range(len(clist)):
        for j in range(i + 1, len(clist)):
            m = NEG
            for u in clist[i]:
                for v in clist[j]:
                    s = sum(map(int.__mul__, u, v))
                    if s > m:
                        m = s
                        if m == LEECH_MAX:
                            break
                else:
                    continue
                break
            cross_max[(i, j)] = cross_max[(j, i)] = m
            cross_attained &= (m == LEECH_MAX)
    rep.check(all(m <= 8 for m in self_max),
              "inside each class the head dot product is at most 8 (cosine 1/4)",
              f"maxima {self_max}")
    rep.check(all(m <= LEECH_MAX for m in cross_max.values()),
              "between different classes it is at most 16 (cosine 1/2)",
              f"maxima {sorted(set(cross_max.values()))}, ceiling {LEECH_MAX} from the "
              f"Leech certification" + ("; every pair attains it" if cross_attained else ""))
    rep.check(all(classes[i].isdisjoint(classes[j])
                  for i in range(len(classes)) for j in range(i + 1, len(classes))),
              f"the {len(classes)} classes are pairwise disjoint")

    # ----------------------------------------------------------- cap-cap
    rep.section("cap-cap pairs")
    worst = None
    ok = True
    for i, di in enumerate(dirs):
        for j, dj in enumerate(dirs):
            tt = rdot(di, dj)
            ci, cj = class_of[di], class_of[dj]
            if i == j:
                m = self_max[ci]              # distinct caps, same direction
            elif ci == cj:
                m = 32                        # the same head on both directions
            else:
                m = cross_max[(ci, cj)]
            lhs = radd((9 * m, 0, 0, 0), tt)
            c = ring_cmp_rational(lhs, Fraction(216))
            if c > 0:
                ok = False
                worst = (i, j, m, rstr(tt))
    # as for the two ring tables below, the threshold must be able to fail:
    # two distinct caps sharing a head sit at 9*32 + t.t', and a pair of
    # directions at cosine +1/2 would give 9*32 + 72 = 360, well above 216.
    # (That is selftest.py's "cap moved to a foreign direction".)
    rep.check(9 * 32 + 72 > 216,
              "the 216 threshold is not vacuous: distinct caps can reach 9*32 + 72",
              f"{9 * 32 + 72} > 216")
    rep.check(ok,
              "9*(head dot) + (tail dot) <= 216 for every pair of directions",
              "" if ok else f"fails at {worst}")
    rep.note("this is exactly <p,q> = (9 y.y' + t.t')/432 <= 1/2 for two caps")
    tails_seen = sorted({rstr(rdot(di, dj)) for di in dirs for dj in dirs})
    rep.note(f"tail dot products occurring: {', '.join(tails_seen)}")

    # ---------------------------------------------------- cap-axis, axis-axis
    rep.section("cap-axis and axis-axis pairs")
    # a cap row has squared norm 432 and an axis row 48, so
    #     <p,q> = t.a / sqrt(432*48) = t.a / 144 <= 1/2  <=>  t.a <= 72.
    rep.check(math.isqrt(432 * 48) ** 2 == 432 * 48,
              "sqrt(432*48) = 144 exactly", "20736 = 144^2")
    THR_CA = Fraction(math.isqrt(432 * 48), 2)
    rep.check(THR_CA == 72, "cap-axis threshold is sqrt(432*48)/2 = 72", str(THR_CA))
    # the threshold must be able to fail: Cauchy-Schwarz allows |t.a| up to
    # sqrt(144*48) = 48 sqrt3 = 83.13..., which exceeds 72, so this is a real test
    rep.check(48 * 48 * 3 > 72 * 72,
              "and that bound is not vacuous: Cauchy-Schwarz permits up to 48 sqrt3 > 72",
              "48^2*3 = 6912 > 5184 = 72^2")
    ok = True
    worst = None
    for d in dirs:
        for a in axis:
            v = rdot(d, a)
            if ring_cmp_rational(v, THR_CA) > 0:
                ok = False
                worst = rstr(v)
    rep.check(ok, "t.a <= 72 for every cap direction and axis point",
              "" if ok else f"fails at {worst}")
    rep.note("this is <p,q> = t.a/144 <= 1/2")
    mx = None
    for d in dirs:
        for a in axis:
            v = rdot(d, a)
            if mx is None or ring_cmp(v, mx) > 0:
                mx = v
    rep.note(f"largest is ({rstr(mx)})/144")

    THR_AA = Fraction(48, 2)
    rep.check(THR_AA == 24, "axis-axis threshold is 48/2 = 24", str(THR_AA))
    rep.check(48 > 24, "and that bound is not vacuous: |a.a'| may reach 48")
    ok = True
    for i, a in enumerate(axis):
        for j, b in enumerate(axis):
            if i == j:
                continue
            v = rdot(a, b)
            if ring_cmp_rational(v, THR_AA) > 0:
                ok = False
    rep.check(ok, "a.a' <= 24 for distinct axis points")
    rep.note("this is <p,q> = a.a'/48 <= 1/2")

    # ------------------------------------------------- equator-* pairs
    rep.section("equator pairs")
    # <p,q> = y.y'/32 <= 1/2  <=>  y.y' <= 16
    rep.check(LEECH_MAX <= 32 // 2,
              "equator-equator: needs y.y' <= 16, and distinct minimal vectors give 16",
              f"bound {LEECH_MAX}, threshold {32 // 2}")
    # <p,q> = 3 y.y'/sqrt(32*432) = y.y'/(16 sqrt6) <= 1/2  <=>  y.y' <= 8 sqrt6
    rep.check(32 * 432 == (48 ** 2) * 6, "sqrt(32*432) = 48 sqrt6", "13824 = 48^2 * 6")
    thr_ec = math.isqrt(8 * 8 * 6)                       # floor(8 sqrt6) = floor(sqrt384)
    rep.check(thr_ec == 19 and 19 * 19 <= 384 < 20 * 20,
              "equator-cap threshold is floor(8 sqrt6) = 19", "361 <= 384 < 400")
    rep.check(LEECH_MAX <= thr_ec,
              "equator-cap: needs y.y' <= 19, and condition (E) makes the heads distinct",
              f"bound {LEECH_MAX}, threshold {thr_ec}")
    rep.note("equator-axis: the inner product is identically 0")

    # ----------------------------------------------------------------- done
    print()
    print(f"  {len(equator)} equator + {len(caps)} caps + {len(axis)} axis = {n} points")
    print(f"  elapsed {time.time() - t0:.1f}s")
    return rep.done(
        f"VERIFIED: {os.path.basename(path)} is a valid kissing configuration of "
        f"{n} points\n           in R^27, so K(27) >= {n}.")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT))
