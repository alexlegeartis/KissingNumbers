#!/usr/bin/env python3
"""
Exact verification of

    tau(38) >= 570236        (previous published record 566652)
    tau(39) >= 756116        (previous published record 755988)

via the Edel-Rains-Sloane construction with support chain (n, 8, 2).

Everything here is integer / exact-rational arithmetic.  No floating point is used
in any decision.

What is verified computationally
--------------------------------
  (A) the two constant-weight codes  A(38,8,8) >= 3025  and  A(39,8,8) >= 3324
      are read from the shipped files and fully checked (size, length, weight,
      minimum distance) by exhaustive pairwise comparison;
  (B) the sign code of level 2 (the even-weight code of length 8, 128 words,
      minimum distance 2) is built explicitly and fully checked;
  (C) the level-3 data (all weight-2 supports, all 4 sign patterns) is explicit;
  (D) every one of the ten geometric conditions of the construction is checked as
      an exact rational inequality, for every pair of levels;
  (E) the arithmetic of the totals.

What is taken from the literature
---------------------------------
  A(38,10) >= 180224 and A(39,10) >= 327680
      (Zinov'ev-Litsyn, "On shortening of codes", Probl. Peredachi Inf. 20:1
      (1984) 3-11; tabulated by Litsyn-Rains-Sloane).
  These are exactly the level-1 codes used by the *existing* records 566652 and
  755988 in H. Cohn's table, so the improvement claimed here does not depend on
  them being correct -- it is the level-2 term that improves, by 128*(3025-2997)
  = 3584 and 128*(3324-3323) = 128 respectively.

Usage:  python scripts/verify.py
"""
import sys, os, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

OK = True
def check(cond, msg):
    global OK
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        OK = False

# ---------------------------------------------------------------- code readers

def read_brouwer(path, n):
    """Brouwer's constant-weight-code format: a '$BASE=16' header then one
    codeword per line as a hex integer whose bit i is coordinate i."""
    words = []
    base = 16
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("$"):
                key, _, val = line[1:].partition("=")
                if key.strip().upper() == "BASE":
                    base = int(val)
                continue
            words.append(int(line, base))
    supports = []
    for w in words:
        s = frozenset(i for i in range(w.bit_length()) if (w >> i) & 1)
        supports.append(s)
    return supports

def check_cw(supports, n, w, d, M, label):
    print(f"[{label}]  constant weight code, expected n={n} w={w} d={d} M={M}")
    check(len(supports) == M, f"file contains exactly {M} codewords (got {len(supports)})")
    check(len(set(supports)) == M, "all codewords distinct")
    check(all(len(s) == w for s in supports), f"every codeword has weight exactly {w}")
    hi = max(max(s) for s in supports)
    check(hi <= n - 1, f"all coordinates lie in 0..{n-1} (highest used = {hi})")
    tmax = 0
    L = list(supports)
    for i in range(len(L)):
        Li = L[i]
        for j in range(i + 1, len(L)):
            t = len(Li & L[j])
            if t > tmax:
                tmax = t
    dmin = 2 * (w - tmax)
    check(tmax <= w // 2, f"max pairwise support intersection = {tmax} <= {w//2}")
    check(dmin >= d, f"minimum Hamming distance = {dmin} >= {d}")
    return tmax

# ------------------------------------------------------------------ sign codes

def even_weight_code(m):
    return [tuple(b) for b in itertools.product((0, 1), repeat=m) if sum(b) % 2 == 0]

def min_distance(words):
    dm = None
    for i in range(len(words)):
        for j in range(i + 1, len(words)):
            d = sum(a != b for a, b in zip(words[i], words[j]))
            if dm is None or d < dm:
                dm = d
    return dm

# ------------------------------------------------------------- the geometry

def geometry(n, c1, c2=8, c3=2):
    """Every level-i vector has support of size c_i and entries +-sqrt(c1/c_i),
    hence squared norm c1.  Two centres at squared distance >= c1 means inner
    product <= c1/2.  Check every case as an exact rational inequality.
    Returns the list of (description, bound_holds)."""
    half = Fraction(c1, 2)
    res = []

    # level 1 vs level 1, same support (the whole c1-set), Hamming distance d1
    d1 = -(-c1 // 4)                      # ceil(c1/4)
    ip = Fraction(c1 - 2 * d1)            # entries are +-1
    res.append((f"L1-L1  d>={d1}: ip = {ip} <= {half}", ip <= half))

    # level 2 vs level 2, same support, sign distance >= ceil(c2/4)
    d2 = -(-c2 // 4)
    ip = Fraction(c1, c2) * (c2 - 2 * d2)
    res.append((f"L2-L2 same support, d>={d2}: ip = {ip} <= {half}", ip <= half))

    # level 2 vs level 2, different supports, |S ^ S'| <= c2//2
    t = c2 // 2
    ip = Fraction(c1, c2) * t
    res.append((f"L2-L2 diff support, |cap|<={t}: ip <= {ip} <= {half}", ip <= half))

    # level 3 vs level 3, same support
    d3 = -(-c3 // 4)
    ip = Fraction(c1, c3) * (c3 - 2 * d3)
    res.append((f"L3-L3 same support, d>={d3}: ip = {ip} <= {half}", ip <= half))

    # level 3 vs level 3, different supports, |cap| <= 1
    ip = Fraction(c1, c3) * 1
    res.append((f"L3-L3 diff support, |cap|<=1: ip <= {ip} <= {half}", ip <= half))

    # cross levels: |ip| <= sqrt(c1/ci)*sqrt(c1/cj)*t = c1*t/sqrt(ci*cj)
    # so the condition c1*t/sqrt(ci*cj) <= c1/2 is  4*t^2 <= ci*cj  (exact, integer)
    for (i, ci), (j, cj) in itertools.combinations([(1, c1), (2, c2), (3, c3)], 2):
        t = min(ci, cj)                   # worst case overlap
        res.append((f"L{i}-L{j} overlap <= {t}: 4*{t}^2 = {4*t*t} <= {ci}*{cj} = {ci*cj}",
                    4 * t * t <= ci * cj))
    return res

# ------------------------------------------------------------------- main

def run(n, c1, A_lvl1, cwfile, M2, tau_expected, tau_old):
    print("=" * 74)
    print(f"DIMENSION {n}:  chain (n0,n1,n2) = ({c1},8,2)")
    print("=" * 74)

    sup = read_brouwer(os.path.join(DATA, cwfile), n)
    check_cw(sup, n, 8, 8, M2, cwfile)

    print("[sign code level 2]  even-weight binary code of length 8")
    E = even_weight_code(8)
    check(len(E) == 128, "128 codewords")
    check(min_distance(E) == 2, "minimum distance 2 = ceil(8/4)")

    print("[sign code level 3]  all 4 sign patterns of length 2")
    T = [tuple(b) for b in itertools.product((0, 1), repeat=2)]
    check(len(T) == 4, "4 codewords")
    check(min_distance(T) == 1, "minimum distance 1 = ceil(2/4)")

    print("[level 3 supports]  all 2-subsets of an n-set")
    n2 = n * (n - 1) // 2
    check(n2 == len(list(itertools.combinations(range(n), 2))), f"C({n},2) = {n2}")

    print("[chain condition]  n >= n0 >= 4*n1 >= 16*n2")
    check(n >= c1 >= 4 * 8 >= 16 * 2, f"{n} >= {c1} >= {4*8} >= {16*2}")

    print("[geometry]")
    for msg, good in geometry(n, c1):
        check(good, msg)

    print("[totals]")
    t1 = A_lvl1 * 1
    t2 = M2 * 128
    t3 = n2 * 4
    tau = t1 + t2 + t3
    print(f"  level 1 : A({c1},{-(-c1//4)}) * 1        = {t1}")
    print(f"  level 2 : A({n},8,8) * A(8,2)  = {M2} * 128 = {t2}")
    print(f"  level 3 : C({n},2)   * A(2,1)  = {n2} * 4 = {t3}")
    check(tau == tau_expected, f"total = {tau} (claimed {tau_expected})")
    print(f"  previous published record : {tau_old}")
    print(f"  improvement               : +{tau - tau_old}")
    print()
    return tau

if __name__ == "__main__":
    run(38, 38, 180224, "a38.8.8.3025H", 3025, 570236, 566652)
    run(39, 39, 327680, "a39.8.8.3324H", 3324, 756116, 755988)
    print("=" * 74)
    print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
    print("=" * 74)
    sys.exit(0 if OK else 1)
