#!/usr/bin/env python3
"""Are the published records in dimensions 9-16 MAXIMAL?

A point can be added to a configuration C iff

    m(C) = min over unit z of max over x in C of <z,x>   <=  1/2,

equivalently iff max |z| over the polytope {z : <z,x> <= 1/2} is at least 1.  Dimension 25
gained +2 exactly this way.  Dimensions 11, 12 and 13 are maximal
(m = 1/sqrt3, 0.604964075, 2/sqrt13); the 0.6049 is for the current 841-point record, not
the 840-point one it replaced.  This sweeps the rest of the low range using Cohn's
own published coordinates, which is the fastest possible check for a free sphere.

Cohn's file lists several configurations per dimension; every one is tested, since a
non-record configuration can still be the one with room (and if a SMALLER configuration turns
out to be extendable by enough points it is worth knowing).

    python sweep.py [dims...]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from extract import blocks, SRC            # noqa: E402
from maximal841 import minmax              # noqa: E402

REC = {9: 306, 10: 510, 11: 604, 12: 841, 13: 1154, 14: 1932, 15: 2564, 16: 4320,
       17: 5730, 18: 7654, 19: 11948, 20: 19448, 21: 29768, 22: 49896, 23: 93150,
       24: 196560}

# m(C) for the RECORD configuration of each dimension this script reaches from Cohn's
# file, to nine decimals.  Dimension 11 is absent because extract.py's grammar does not
# parse that block; dimensions 17, 18 and 19 come from separate runs (see README.md).
EXPECT = {9: 0.577350269, 10: 0.670820393, 12: 0.604964075, 13: 0.554700196,
          14: 0.534522484, 15: 0.534522484, 16: 0.577350269}
TOL = 1e-6

if __name__ == '__main__':
    # a .txt argument is the path to Cohn's data set, consumed by extract._find_source
    want = [int(a) for a in sys.argv[1:] if not a.endswith('.txt')] or list(range(9, 17))
    print(__doc__)
    got = {}
    try:
        for b in blocks(SRC):
            if b['dim'] in want:
                got.setdefault(b['dim'], []).append(b)
    except SystemExit as e:
        print(e)
        print("SKIP -- Cohn's data set is not present; nothing checked.")
        raise SystemExit(0)
    seen = {}
    for dim in sorted(got):
        print("=" * 70)
        for b in got[dim]:
            X = np.array(b['pts'])
            if len(X) != b['n']:
                print("dim %2d: parsed %d of %d -- SKIP" % (dim, len(X), b['n']))
                continue
            c = np.array(b['c'])
            Y = X * np.sqrt(c)[None, :]
            U = Y / np.sqrt((Y * Y).sum(1))[:, None]
            G = U @ U.T
            np.fill_diagonal(G, -1.0)
            tag = "RECORD" if b['n'] == REC.get(dim) else "alternate"
            print("dim %2d, %d points (%s): max inner product %.9f"
                  % (dim, b['n'], tag, G.max()))
            if G.max() > 0.5 + 1e-7:
                print("    not a valid 60-degree code at this tolerance -- SKIP")
                continue
            best, bz = minmax(U, restarts=120, seed=5)
            m = 0.5 / best if best > 0 else 9
            print("    max |z| = %.9f -> m(C) = %.9f  %s"
                  % (best, m, "FREE POINT!" if best >= 1 - 1e-9 else
                     "maximal (margin %.1f%%)" % (100 * (m - 0.5) / 0.5)))
            if tag == "RECORD":
                seen[dim] = min(seen.get(dim, 9.0), m)
            if best >= 1 - 1e-9:
                z = bz / np.linalg.norm(bz)
                print("    *** tau(%d) >= %d   max <z,x> = %.12f"
                      % (dim, b['n'] + 1, (U @ z).max()))
                np.save(os.path.join(HERE, 'free_%d_%d.npy' % (dim, b['n'])), z)

    # ---- verdict.  The optimisation is numerical; this pins its output to the published
    # ---- values, so a regression in the ascent or in the parser cannot pass unnoticed.
    print("=" * 70)
    bad = []
    for dim in sorted(d for d in EXPECT if d in want):
        if dim not in seen:
            bad.append("dim %d: the record configuration was never reached" % dim)
            continue
        m = seen[dim]
        if m <= 0.5:
            bad.append("dim %d: m(C) = %.9f <= 1/2, so the record is NOT maximal" % (dim, m))
        elif abs(m - EXPECT[dim]) > TOL:
            bad.append("dim %d: m(C) = %.9f, expected %.9f" % (dim, m, EXPECT[dim]))
        else:
            print("  dim %2d  m(C) = %.9f  maximal, margin %.1f%%"
                  % (dim, m, 100 * (m - 0.5) / 0.5))
    print("=" * 70)
    if bad:
        for line in bad:
            print("  FAIL  " + line)
        print("*** THE SWEEP DISAGREES WITH THE RECORDED VALUES ***")
        raise SystemExit(1)
    print("every record configuration reached is maximal, at the recorded m(C)")
