"""Independent verification of the dimension-19 reduction, from the graph itself.

The claim chain is:
    tau(19) = 10668 + 4 * alpha(Cay(K,S))            [section 31, four isomorphic components]
    alpha(Cay(K,S)) = maxbip(G_0) <= 2 * alpha(G_0)  [section 71]
    alpha(G_0) = 160, proven optimal                 [alpha512.py]
=>  alpha(Cay(K,S)) <= 320, achieved by Ho  =>  tau(19) = 11948 is the exact ceiling.

The one structural step that is not a solver output is the split, so it is checked here
directly on the 1024-vertex graph rather than argued: that S has exactly one element outside
the even-weight hyperplane K_0, that each half induces G_0, and that the cross edges are a
PERFECT MATCHING.  Everything is checked by explicit enumeration.
"""
import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
import numpy as np
from capcheck import golay

GOL = np.array(sorted(golay()), dtype=np.uint8)[:, 5:]
W = GOL.sum(1)
codes = [int("".join(map(str, r)), 2) for r in GOL]
Sraw = sorted(codes[i] for i in range(len(codes)) if 0 < W[i] <= 4)
print(__doc__)
print("S has %d elements, weights %s" % (len(Sraw), dict(collections.Counter(
    bin(c).count('1') for c in Sraw))))

span = {0}
for s in Sraw:
    span |= {x ^ s for x in span}
K = sorted(span)
print("K = span(S) has %d elements (dimension %d)" % (len(K), len(K).bit_length() - 1))
Sset = set(Sraw)

K0 = [v for v in K if bin(v).count('1') % 2 == 0]
odd = [s for s in Sraw if bin(s).count('1') % 2 == 1]
print("K_0 (even-weight part of K): %d elements" % len(K0))
assert len(K0) * 2 == len(K), "K_0 is not a hyperplane of K"
assert len(odd) == 1, "S has %d odd-weight elements, expected 1" % len(odd)
t = odd[0]
assert t not in set(K0)
S0 = [s for s in Sraw if s != t]
assert all(s in set(K0) for s in S0), "some weight-4 word is not in K_0"
print("t = %d is the unique element of S outside K_0 ; |S_0| = %d" % (t, len(S0)))

K0set = set(K0)
half2 = {v ^ t for v in K0}
assert not (K0set & half2) and len(K0set | half2) == len(K), "K is not K_0 u (K_0+t)"

# each half induces G_0
for v in K0:
    nb = {v ^ s for s in Sraw}
    inside = nb & K0set
    assert inside == {v ^ s for s in S0}, "half 1 does not induce G_0 at %d" % v
print("both halves induce Cay(K_0, S_0): verified on all %d vertices" % len(K0))

# cross edges are a perfect matching
bad = 0
for v in K0:
    cross = {u for u in ({v ^ s for s in Sraw}) if u in half2}
    if cross != {v ^ t}:
        bad += 1
assert bad == 0, "%d vertices have cross-degree != 1" % bad
print("cross edges form a PERFECT MATCHING x <-> x^t: verified on all %d vertices" % len(K0))

deg = {len({v ^ s for s in S0}) for v in K0}
print("G_0 is %s-regular on %d vertices" % (deg, len(K0)))
print("\nreduction verified:  alpha(Cay(K,S)) = maxbip(G_0) <= 2*alpha(G_0)")
print("with alpha(G_0) = 160 (CP-SAT, OPTIMAL):  alpha(Cay(K,S)) <= 320")
print("tau(19) <= 10668 + 4*320 = %d, and Ho attains it => tau(19) = 11948 is the CEILING"
      % (10668 + 4 * 320))
