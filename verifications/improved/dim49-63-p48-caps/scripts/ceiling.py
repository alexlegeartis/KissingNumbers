#!/usr/bin/env python3
"""The absolute ceiling of this construction, dimension by dimension.

    tau(48+k) >= 52 416 000 + 2*S(lambda(k)) + tau(k)

Everything on the right except S is fixed: 52 416 000 is tau(P_48), tau(k) is the pole layer,
and lambda(k) is half the largest antipodally symmetric 60-degree code of R^k, which is capped
by the kissing number record in dimension k.  So the only free quantity is S, the total size of
the lambda(k) disjoint classes, and S is capped by lambda(k) times the largest class that
exists.  That gives a hard ceiling

    52 416 000 + 2*lambda(k)*Cmax + tau(k),      Cmax = 7069,

which no improvement to the class family can exceed.  Cmax = 7069 is itself a ceiling, not just
the best found: it is the independence number of the conflict graph on the P_48 minimal lines
(see ../dim49-63-p48-caps/README.md and KNOWLEDGE 89).

WHY THIS MATTERS.  A single sign code is already a kissing configuration:

    tau(n) >= max over n_0 <= n of A(n_0, ceil(n_0/4)),

and the codes [62,26,16] and [63,27,16] give 2^26 = 67 108 864 and 2^27 = 134 217 728, with
ceil(62/4) = ceil(63/4) = 16.  Comparing those against the ceilings below settles whether this
construction can ever compete in those two dimensions.  It cannot.
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CMAX = 7069
TAU_P48 = 52416000
LINES = 26208000        # minimal LINES of P_48, half of tau(P_48) = 52416000
# largest k with minimum distance >= ceil(n_0/4) among best known binary linear codes
# (codetables.de, read 2026-08-23); a group of four carries down by shortening
ERS0 = {48: 24, 52: 21, 56: 23, 60: 25, 61: 25, 62: 26, 63: 27}


def main():
    rows = json.load(open(os.path.join(HERE, '..', 'rows49-63.json')))
    sizes = np.load(os.path.join(HERE, '..', 'data', 'class_sizes_p48.npy'))
    assert int(sizes.max()) == CMAX, (int(sizes.max()), CMAX)
    print(__doc__)
    print("  largest class in the shipped family: %d  (the ceiling Cmax)" % int(sizes.max()))
    print()
    print("  dim   k  lambda   claimed      ceiling      ERS level 0   verdict")
    dead = []
    for r in rows:
        d, k, lam = r['dim'], r['k'], r['lam']
        tau_k = r['gain'] - 2 * r['S']          # gain = 2S + tau(k), the pole layer
        ceil_ = TAU_P48 + 2 * lam * CMAX + tau_k
        assert r['value'] <= ceil_, (d, r['value'], ceil_)
        lvl0 = max([1 << kk for n0, kk in ERS0.items() if n0 <= d] or [0])
        if ceil_ <= lvl0:
            verdict = 'DEAD -- ceiling below the level-0 floor'
            dead.append((d, ceil_, lvl0))
        elif lvl0 >= r['value']:
            verdict = 'claim beaten, but the ceiling is above it'
        else:
            verdict = 'clear'
        print("   %2d  %2d  %5d   %-11d  %-11d  %-11d   %s"
              % (d, k, lam, r['value'], ceil_, lvl0, verdict))
    print()
    print("  HOW MUCH OF THE GAP TO THE CEILING IS REACHABLE")
    print("  Classes are random automorphism images of one base class, and two images of a")
    print("  class of size C share about C^2/L lines by chance, L = %d being the number of"
          % LINES)
    print("  minimal LINES.  Summing that over the lambda images is the overlap one must expect")
    print("  from randomness alone; best-of-N selection attacks only the variance around it,")
    print("  worth about 2.1 standard deviations per class.")
    print()
    print("  The gap is quoted in LINES here, not in the bound, since the bound counts each")
    print("  line twice.")
    print()
    print("  dim   lambda   gap in lines   expected overlap   % of gap   best-of-N can reach")
    for r in rows:
        d, lam = r['dim'], r['lam']
        tau_k = r['gain'] - 2 * r['S']
        gap = ((TAU_P48 + 2 * lam * CMAX + tau_k) - r['value']) // 2   # bound counts twice
        exp = sum(CMAX * CMAX * k / float(LINES) for k in range(lam))
        var = sum(2.1 * (CMAX * CMAX * k / float(LINES)) ** 0.5 for k in range(lam))
        if gap <= 0:
            continue
        print("   %2d   %5d   %11d   %16.0f   %7.1f%%   %19.0f"
              % (d, lam, gap, exp, 100 * exp / gap, var))
    print()
    print("  A figure above 100% is not an error: it means the family is already doing better")
    print("  than plain random images, which is the best-of-N selection in regenerate.py")
    print("  earning its keep.  Everywhere the figure is near 100%, the distance to the ceiling")
    print("  is explained by chance overlap and nothing else.")
    print()
    print("  So the family is within a few tens of thousands of lines of what random images can")
    print("  give.  Closing the rest needs a family disjoint BY CONSTRUCTION -- a subgroup that")
    print("  acts freely on the base class, or a coset decomposition -- not more candidates.")
    print()
    if dead:
        print("  Dimensions %s cannot be rescued by ANY class family: even with every one of"
              % ', '.join(str(d) for d, _, _ in dead))
        print("  the lambda(k) classes at the maximum size 7069 the total stays below what a")
        print("  single sign code already gives.  They are not improvements and cannot become")
        print("  improvements; raising lambda(k) would itself require a new kissing record in")
        print("  dimension k.")
        for d, c, l in dead:
            print("     dim %d: ceiling %d < %d, short by %d" % (d, c, l, l - c))
    else:
        print("  every dimension's ceiling is above its level-0 floor")


if __name__ == '__main__':
    main()
