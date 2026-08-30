#!/usr/bin/env python3
"""Stamp each record axis layer with a digest of the source it indexes into.

    python seal_record.py

The shipped file carries an index set `keep` into the source configuration, and the source is
REBUILT by the verifier rather than stored -- from the cap file or the Golay code, the odd
sign convention, and the shipped code words.  That is the right way round, but it makes the
ORDER of the rebuilt rows part of the contract: if it ever differed from the order the search
used, `keep` would select different poles and the verification would be checking something
else.  A wrong subset would almost certainly fail the cap test, but "almost certainly" is not
what the rest of this package settles for.

So the digest of the source, as the search built it, goes in the file, and
scripts/verify_record.py recomputes it from its own rebuild.  Running this file also
re-derives every layer from scratch and checks that the shipped index set is exactly the set
of poles that clear -- so it is a verification in its own right, not only a stamp.
"""
import hashlib
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from axis_record import cleared, source                            # noqa: E402
from axis_rotate import directions, rows                           # noqa: E402
import axis_cohn                                                   # noqa: E402

DATA = os.path.join(HERE, '..', 'data')


def digest(T):
    """sha256 of the source, in a canonical layout"""
    A = np.ascontiguousarray(T, dtype=np.int64)
    return np.frombuffer(hashlib.sha256(A.tobytes()).digest(), dtype=np.uint8)


def main():
    print(__doc__)
    n = 0
    for k in range(9, 24):
        p = os.path.join(DATA, 'poles_rec_%d.npz' % k)
        if not os.path.exists(p):
            continue
        z = dict(np.load(p))
        M = z['M'].astype(np.int64)
        if 'cs' in z:
            # a PUBLISHED source, in its own metric: scripts/axis_cohn.py.  The source is
            # read rather than rebuilt, so what this re-derives is the index set.
            cs = z['cs'].astype(np.int64)
            T, c = axis_cohn.source(k)
            assert (c == cs).all(), 'k=%d: the shipped source has a different metric' % k
            W, _cw, m = directions(k, rows())
            mu = int((M.T @ M)[0, 0]) // int(cs[0])
            good, _w = axis_cohn.cleared(T, cs, z['N'].astype(np.int64), int(z['D'][0]),
                                         W @ M, mu, m)
        else:
            code = [int(c) for c in z['code']]
            Zk, T, split, code2, M2, mu, oct_, ncol = source(k, code=code, M=M)
            assert code2 == sorted(code)
            good = cleared(T, Zk, z['N'].astype(object), int(z['D'][0]))
        idx = np.nonzero(good)[0].astype(np.int32)
        assert len(idx) == len(z['keep']), \
            'k=%d: the layer re-derives as %d poles, the file says %d' % (k, len(idx),
                                                                         len(z['keep']))
        assert (idx == z['keep']).all(), 'k=%d: the index set differs' % k
        z['sha'] = digest(T)
        np.savez(p, **z)
        print('   k = %2d: %5d poles re-derived and matched; source digest %s'
              % (k, len(idx), ''.join('%02x' % b for b in z['sha'][:8])))
        n += 1
    print()
    print('   sealed %d layer(s)' % n)
    return 0


if __name__ == '__main__':
    sys.exit(main())
