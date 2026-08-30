"""alpha(Cay(K', W_4)) = 192 EXACTLY, so tau(17) = 5730 in the layered family.

    python theta17.py

This closes what sections 74-75 of KNOWLEDGE.md left open.  It replaces the CP-SAT runs
entirely: those were a lottery (three runs of the same capped model returned dual bounds 192,
193 and 195, and the deterministic single-worker variant returned 238).  Nothing here is a
search.  The whole proof is one linear-programming DUAL, exhibited as six rational numbers,
and checking it is integer arithmetic on 1024 values.

WHY AN LP DECIDES IT.  The obstruction all along was that K' was being treated as an arbitrary
1024-vertex graph.  It is not: K' is an elementary abelian group -- a linear code, closed under
XOR -- and W_4 is a symmetric subset of it, so G = Cay(K', W_4) is a translation scheme and
Delsarte's linear programming bound applies verbatim.  For a Cayley graph on an abelian group
that bound equals Lovasz's theta', so it is far stronger than anything the fibre relaxation can
say (the clique cover forces 256, and so does the Hoffman ratio bound, since lambda_min = -20
gives exactly 1024*20/80 = 256).  The LP gives 192, on the nose.

THE CERTIFICATE.  Let C be an independent set, so no difference of two of its elements lies in
W_4.  For any g : K' -> Q whose Fourier transform is non-negative,

    sum_{c,c' in C} g(c-c')  =  (1/1024) sum_y ghat(y) |sum_{c in C} chi_y(c)|^2  >=  0 ,

and if also g(x) <= -1 for every x outside W_4 u {0} the left side is at most
|C|*g(0) - |C|(|C|-1).  Hence |C| <= 1 + g(0).  We exhibit such a g with g(0) = 191.

It is given by its Fourier coefficients, which take only six values.  Index the characters y of
K' by two invariants: the eigenvalue lambda(y) = sum_{w in W_4} chi_y(w), and whether y is
trivial on the fibre subgroup H (the fibre through 0).  Then

    b_y = ghat(y)/1024 =    65/80  on the  16 characters with (lambda,triv) = ( 30, no)
                            40/80  on the  10             with              = ( 12, yes)
                            30/80  on the 120             with              = ( 12, no)
                            28/80  on the 160             with              = (  6, no)
                            10/80  on the 360             with              = ( -4, no)
                             9/80  on the 240             with              = ( -2, no)
                             0     on the remaining 118

and the masses come to 13 + 5 + 45 + 56 + 45 + 27 = 191 exactly.  Everything below is checked
in integers scaled by 80; no floating point is used and no solver is called.

THE LOWER BOUND is Cohn-Li's own configuration: six pairwise non-adjacent fibres, 192 points,
every one of the 18336 pairs checked.  So alpha = 192 and

    tau(17) = 5346 + 2*192 = 5730 .

SCOPE.  Unchanged from section 75: this is the layered family over a Lambda_16 equator with an
antipodal pole pair and integral tiers -- the family containing Cohn-Li's construction and every
usable variant of it (section 74i).  It is not a proof that tau(17) <= 5730 for arbitrary
configurations.  What it does settle is that no configuration IN THAT FAMILY beats 5730, which
is what the repository listed as its one live lead.
"""
import collections
import itertools
import os
import sys
from fractions import Fraction as Fr

import numpy as np

DEN = 80
CERT = {(30, False): 65, (12, True): 40, (12, False): 30,
        (6, False): 28, (-4, False): 10, (-2, False): 9}

H = os.path.dirname(os.path.abspath(__file__))


def build():
    """K', W_4 and the fibres, rebuilt from the supports exactly as cap22.py does."""
    sup = [tuple(int(v) for v in r)
           for r in np.load(os.path.join(H, '..', 'data', 'd17_supports.npy'))]
    rm1 = {0}
    for f in range(1, 16):
        for c in (0, 1):
            v = sum(1 << x for x in range(16) if bin(x & f).count('1') % 2 == c)
            rm1 |= {y ^ v for y in rm1}
    rm2 = [x for x in range(1 << 16) if all(bin(x & y).count('1') % 2 == 0 for y in rm1)]
    m0 = sum(1 << i for i in sup[0])
    kl = sorted(b for b in rm2 if bin(b & m0).count('1') % 2 == 0)
    w4 = [b for b in kl if bin(b).count('1') == 4]
    fib, reps = {}, []
    for b in kl:
        if b not in fib:
            for r in rm1:
                fib[b ^ r] = len(reps)
            reps.append(b)
    return kl, w4, fib, len(reps)


def coordinates(kl):
    """K' is closed under XOR; return an F_2-linear coordinate map onto F_2^10."""
    red = {}
    for b in kl:
        v = int(b)
        for p in sorted(red, reverse=True):
            if v >> p & 1:
                v ^= red[p]
        if v:
            red[v.bit_length() - 1] = v
    piv = sorted(red)

    def co(b):
        v, out = int(b), 0
        for k, p in enumerate(piv):
            if v >> p & 1:
                v ^= red[p]
                out |= 1 << k
        assert v == 0, "K' is not spanned by the pivots"
        return out

    return co, len(piv)


def wht(vec):
    """Walsh-Hadamard transform, exact: returns W(x) = sum_y vec[y] (-1)^<x,y>."""
    a = list(vec)
    n = len(a)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                a[j], a[j + h] = a[j] + a[j + h], a[j] - a[j + h]
        h *= 2
    return a


def wht_bruteforce(vec):
    """The same transform by definition, as a control on the butterfly above.

    A fast Walsh-Hadamard is easy to get subtly wrong -- a swapped sign or a mis-strided
    butterfly still returns something that looks like a spectrum.  The whole bound rests on
    g, so g is recomputed here straight from the character sums.  1024 x 1024 integer terms,
    about a second, and it removes the transform from the trusted base entirely.
    """
    n = len(vec)
    return [sum(vec[y] * (1 if bin(x & y).count("1") % 2 == 0 else -1) for y in range(n))
            for x in range(n)]


def certificate(kl, w4, hsub, show=None):
    """alpha(Cay(kl, w4)) <= 1 + g(0), checked in integers.  Returns (bound, spectrum).

    Split out so that frame17.py can apply the SAME six numbers to the graph it measures from
    Cohn's coordinates.  That matters: the two scripts label the 16 coordinates differently --
    theta17 builds RM(2,4) in the Boolean-cube order, Cohn's file is in its own order -- so the
    two K's are isomorphic but not equal as sets of words.  Certifying both settles it without
    anyone having to find the permutation.
    """
    co, dim = coordinates(kl)
    assert dim == 10
    d = 1 << dim
    wc = sorted(co(w) for w in w4)
    hs = sorted(co(x) for x in hsub)
    assert len(wc) == 60 and 0 not in wc and len(hs) == 32
    wset = set(wc)
    lam = {y: sum(1 if bin(y & w).count('1') % 2 == 0 else -1 for w in wc) for y in range(d)}
    triv = {y: all(bin(y & h).count('1') % 2 == 0 for h in hs) for y in range(d)}
    b = [CERT.get((lam[y], triv[y]), 0) for y in range(d)]
    assert all(v >= 0 for v in b), "the Fourier coefficients must be non-negative"
    g = wht(b)
    assert g == wht_bruteforce(b), 'the fast transform disagrees with the character sums'
    assert sum(v * v for v in g) == d * sum(v * v for v in b), 'Parseval fails'
    off = [x for x in range(d) if x != 0 and x not in wset]
    assert len(off) == 963
    worst = max(g[x] for x in off)
    if show is not None:
        hist = collections.Counter(g[x] for x in off)
        print("%s: b >= 0 on all %d characters (%d carry mass); g outside W_4 u {0}:"
              % (show, d, sum(1 for v in b if v)))
        for v in sorted(hist, reverse=True):
            print("      g = %-8s at %4d points" % (Fr(v, DEN), hist[v]))
        print("   max = %s  <= -1 : %s" % (Fr(worst, DEN), worst <= -DEN))
    assert worst <= -DEN, "certificate infeasible"
    assert g[0] == sum(b)
    return 1 + Fr(g[0], DEN), collections.Counter(lam.values())


def main():
    print(__doc__)
    kl, w4, fib, nfib = build()
    print("K' : %d elements, W_4 : %d, fibres : %d" % (len(kl), len(w4), nfib))
    assert (len(kl), len(w4), nfib) == (1024, 60, 32)

    # K' really is a group under XOR -- the whole argument rests on this
    ks = set(kl)
    assert all((a ^ b) in ks for a in kl for b in w4), "K' is not closed under XOR"
    assert all((0 ^ w) in ks and w in ks for w in w4)
    print("K' is closed under XOR (checked on all %d x %d pairs) -- a translation scheme"
          % (len(kl), len(w4)))

    co, dim = coordinates(kl)
    assert dim == 10
    d = 1 << dim
    wc = sorted(co(w) for w in w4)
    wset = set(wc)
    assert len(wset) == 60 and 0 not in wset
    hsub = sorted(co(b) for b in kl if fib[b] == fib[kl[0]])
    assert len(hsub) == 32

    # --- the spectrum, for the record: the ratio bound is exactly the clique-cover 256 -------
    lam = {y: sum(1 if bin(y & w).count('1') % 2 == 0 else -1 for w in wc) for y in range(d)}
    spec = collections.Counter(lam.values())
    print("\nspectrum of G  " + "  ".join("%d^%d" % (l, spec[l]) for l in sorted(spec, reverse=True)))
    lmin, deg = min(spec), max(spec)
    print("Hoffman ratio bound : n(-l)/(d-l) = %d*%d/%d = %s  -- no better than the clique cover"
          % (d, -lmin, deg - lmin, Fr(d * -lmin, deg - lmin)))

    # --- the certificate ----------------------------------------------------------------------
    print("")
    hfib = [x for x in kl if fib[x] == fib[kl[0]]]
    bound, _ = certificate(kl, w4, hfib, show="certificate")
    print("   g(0) = %s, so alpha <= 1 + g(0) = %s" % (bound - 1, bound))
    assert bound == 192

    # --- the matching construction ------------------------------------------------------------
    w4s = set(w4)
    qadj = [[0] * nfib for _ in range(nfib)]
    for a in kl:
        for w in w4:
            c, e = fib[a], fib[a ^ w]
            if c != e:
                qadj[c][e] = qadj[e][c] = 1
    six = next(s for s in itertools.combinations(range(nfib), 6)
               if all(qadj[p][q] == 0 for p, q in itertools.combinations(s, 2)))
    pts = [a for a in kl if fib[a] in six]
    assert len(pts) == 192
    bad = sum(1 for a, e in itertools.combinations(pts, 2) if (a ^ e) in w4s)
    print("\nsix pairwise non-adjacent fibres %s -> %d points, %d adjacent pairs among the %d"
          % (list(six), len(pts), bad, len(pts) * 191 // 2))
    assert bad == 0
    seven = any(all(qadj[p][q] == 0 for p, q in itertools.combinations(s, 2))
                for s in itertools.combinations(range(nfib), 7))
    assert not seven

    print("\n" + "=" * 78)
    print("alpha(Cay(K', W_4)) = 192 exactly    (192 <= alpha <= 192)")
    print("tau(17) = 5346 + 2*192 = 5730        -- the layered family is closed")
    print("=" * 78)
    return 0


if __name__ == '__main__':
    sys.exit(main())
