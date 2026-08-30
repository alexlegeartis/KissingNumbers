#!/usr/bin/env python3
"""Independent check that the pool a certificate was built from has no repeated LINE.

    python checkpool.py cert_prod.json

Why this is needed.  verify_classes proves the classes are disjoint as sets of POOL INDICES.
That is disjointness as sets of lines only if the pool itself has no line twice, and the pool
is deduplicated by a 64-bit hash -- so the guarantee is one-sided by construction (a hash
collision DISCARDS a line, it can never let a duplicate through, since equal rows always hash
equal).  This script rebuilds the pool from the seed and re-hashes every row with a FRESH,
independently drawn pair of hash vectors, then asserts all keys are distinct.  An error in the
deduplicator would have to survive two independent hashes to escape.

It also re-checks that every row is a minimal vector and sign-normalised, so that no row is
the negative of another.
"""
import json, os, sys
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')
import numpy as np
import torch
import gamma72 as g7

if len(sys.argv) < 2:
    raise SystemExit("usage: python checkpool.py <cert.json>\n                     independently re-hash the pool and prove it holds no line twice")
cert = json.load(open(sys.argv[1]))
a = cert['args']
g7.select_backend(g7.device())
P = g7.spread_pool(a['pool'], walk=a['walk'], seed=a['seed'], orbit=a['orbit'],
                   nu=a['nu'], sample=a['sample'], maxrows=a.get('maxrows', 24_000_000),
                   verbose=False)
n = len(P)
assert n == cert['pool'], "pool size %d != certified %d" % (n, cert['pool'])
print("pool rebuilt from seed %d: %d lines" % (a['seed'], n), flush=True)

rng = np.random.default_rng(987654321)                    # nothing to do with the run's seed
r1 = rng.integers(1, 2 ** 40, size=72)
r2 = rng.integers(1, 2 ** 40, size=72)
keys = torch.empty(n, dtype=torch.int64, device=P.device)
for a0 in range(0, n, 8_000_000):
    keys[a0:a0 + 8_000_000] = g7._hash_keys(P[a0:a0 + 8_000_000], r1, r2)
u = torch.unique(keys)
print("independent 64-bit re-hash: %d distinct keys out of %d rows" % (len(u), n), flush=True)
assert len(u) == n, "the pool contains a repeated line"

G = torch.tensor(g7.load_gram(), dtype=torch.float64, device=P.device)
for a0 in range(0, n, 2_000_000):
    B = P[a0:a0 + 2_000_000].to(torch.float64)
    assert bool((((B @ G) * B).sum(1) == 8).all()), "row of norm != 8"
    nz = B != 0
    first = torch.argmax(nz.to(torch.uint8), dim=1)
    lead = B.gather(1, first[:, None]).squeeze(1)
    assert bool((lead > 0).all()), "row is not sign-normalised, so +-x could both be present"
print("every row minimal (norm 8) and sign-normalised")
print("POOL IS DUPLICATE-FREE -- disjointness by pool index is disjointness as lines")

# The pool this script just rebuilt is only THE pool if its fingerprint matches the one the
# certificate was written against; see gamma72.pool_digest.
d = g7.pool_digest(P)
if cert.get('pool_digest'):
    ok = d == cert['pool_digest']
    print("pool digest %s  certificate %s  ->  %s"
          % (d[:16], cert['pool_digest'][:16], "MATCH" if ok else "*** DIFFERENT POOL ***"))
    assert ok, "this is not the pool the certificate was written against"
else:
    print("pool digest %s  (the certificate predates pool digests and carries none)" % d[:32])
