#!/usr/bin/env python3
"""Zero-sum triple partition of the 196560 minimal vectors of the Leech lattice.

Vectors are in Cohn's scaling, norm^2 = 32, so inner products lie in {0,+-8,+-16,+-32}
and cos = <a,b>/32.  A zero-sum triple is {a,b,c} with a+b+c = 0; then <a,b> = -16,
i.e. cos = -1/2, exactly the 120-degree condition the cap construction needs at t = 2/3.

Greedy in index order.  Output: triples24.npy of shape (T,3), indices into leech196560.npy.
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from leech import load

A = load().astype(np.int32)
n = len(A)
print("Leech minimal vectors:", A.shape, "norms", set((A**2).sum(1).tolist()))

# exact hash: vectors have entries in [-4,4]; pack base 9
POW = (9 ** np.arange(24, dtype=np.int64))
key = ((A.astype(np.int64) + 4) * POW).sum(1)
assert len(np.unique(key)) == n
order = np.argsort(key)
skey = key[order]

def index_of(vecs):
    """row indices of vecs in A, or -1"""
    k = ((vecs.astype(np.int64) + 4) * POW).sum(1)
    p = np.searchsorted(skey, k)
    p = np.clip(p, 0, n - 1)
    hit = skey[p] == k
    out = np.where(hit, order[p], -1)
    return out

Af = A.astype(np.float32)
used = np.zeros(n, dtype=bool)
triples = []
t0 = time.time()
for i in range(n):
    if used[i]:
        continue
    d = np.rint(Af @ Af[i]).astype(np.int32)
    cand = np.flatnonzero((d == -16) & ~used)
    if len(cand) == 0:
        continue
    C = -A[i] - A[cand]
    m = index_of(C)
    ok = (m >= 0) & ~used[np.maximum(m, 0)] & (m != i) & (m != cand)
    w = np.flatnonzero(ok)
    if len(w) == 0:
        continue
    j = int(cand[w[0]]); k3 = int(m[w[0]])
    used[i] = used[j] = used[k3] = True
    triples.append((i, j, k3))
    if len(triples) % 5000 == 0:
        print("  %6d triples, %6d used, %.0fs" % (len(triples), used.sum(), time.time() - t0), flush=True)

T = np.array(triples, dtype=np.int32)
print("triples:", len(T), " covered:", 3 * len(T), "of", n, " leftover:", n - 3 * len(T))
np.save(os.path.join(HERE, 'data', 'triples24.npy'), T)
print("saved triples24.npy")
