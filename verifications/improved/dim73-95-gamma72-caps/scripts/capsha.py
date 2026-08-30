#!/usr/bin/env python3
"""The identity of a cap direction set, as a digest, so a layer cannot be read against
the wrong one.

    python capsha.py        print the digest of every shipped cap set

WHY.  An axis layer is stored as a rotation N/D together with a set of INDICES into the cap
direction set W.  Both halves are meaningless against a different W: the indices point at
vectors that may not be there, and a rotation that clears one set of caps need not clear
another.  Until 2026-08-30 that could not go wrong, because each k had exactly one cap set
ever chosen for it.  Then k = 12 and 13 acquired a second (the Kappa sections K_12 and K_13,
which are larger than Lambda_12 and Lambda_13 and split perfectly into zero-sum triples), and
the old layers stayed on disk indexing vectors that the dimension no longer uses.

Nothing numeric catches that.  A layer of 648 poles read against a 756-vector cap set is a
perfectly well-typed thing to compute; it is simply an answer to a question nobody asked.  So
a layer written from now on carries the digest of the W it was built from, and a reader that
finds a digest checks it.  A layer without one is from before this file and is accepted, which
is why the file names were kept apart as well: poles_kap_<k>.npz for the Kappa cap sets,
poles_rot_<k>.npz for the others.

The digest covers the vectors AND the metric, since <x,y> = sum_i c_i x_i y_i and a set of
coordinates without its c is not a configuration.  Row order is part of it: the indices in a
layer are positions in W, so a reordering of W is a different object for this purpose even
though it is the same configuration.
"""
import hashlib
import os
import sys

import numpy as np


def digest(W, c):
    """the identity of a cap set: its vectors, in order, and its metric"""
    Wi = np.ascontiguousarray(np.asarray(W, dtype=np.int64))
    ci = np.ascontiguousarray(np.asarray(c, dtype=np.int64))
    assert Wi.dtype == np.int64 and ci.dtype == np.int64, 'the digest needs plain integers'
    h = hashlib.sha256()
    h.update(b'capset\x00')
    h.update(np.array(Wi.shape, dtype=np.int64).tobytes())
    h.update(Wi.tobytes())
    h.update(ci.tobytes())
    return np.frombuffer(h.digest(), dtype=np.uint8).copy()


def matches(z, W, c):
    """True if the layer file z was built from this cap set, or carries no claim at all"""
    if 'wsha' not in getattr(z, 'files', []):
        return True
    return bool((np.asarray(z['wsha'], dtype=np.uint8) == digest(W, c)).all())


if __name__ == '__main__':
    DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
    print(__doc__)
    print('   cap set              |W|  x dim   digest (first 16 hex)')
    n = 0
    for k in range(9, 24):
        for wf, cf in (('lam%d_W.npy' % k, None), ('kap%d_W.npy' % k, 'kap%d_c.npy' % k),
                       ('rec_%d_W.npy' % k, 'rec_%d_c.npy' % k),
                       ('lambda9_W.npy' if k == 9 else '', None)):
            if not wf:
                continue
            p = os.path.join(DATA, wf)
            if not os.path.exists(p):
                continue
            W = np.load(p).astype(np.int64)
            c = (np.load(os.path.join(DATA, cf)).astype(np.int64) if cf
                 else np.ones(W.shape[1], np.int64))
            d = digest(W, c)
            print('   %-20s %6d x %-4d %s' % (wf, W.shape[0], W.shape[1],
                                              ''.join('%02x' % v for v in d[:8])))
            n += 1
    print()
    print('   %d cap sets' % n)
