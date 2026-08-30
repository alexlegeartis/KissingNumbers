#!/usr/bin/env python3
"""Rebuild the classes from a certificate and check they match, bit for bit.

    python regenerate.py cert_prod.json

Reads only the seed and the parameters, reruns the identical pipeline, re-verifies the classes
in exact arithmetic and compares the SHA-256.  Runs on CPU if there is no GPU -- much slower,
but the result must be identical, which is the whole point of taking every random decision in
numpy's PCG64 and leaving the device nothing but exact integer GEMM.
"""
import hashlib, json, os, sys
os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'expandable_segments:True')
import numpy as np
import gamma72 as g7

if len(sys.argv) < 2:
    raise SystemExit("usage: python regenerate.py <cert.json> [n]\n"
                     "       rebuild the classes from the certificate's seed and match its "
                     "SHA-256;\n       an optional n checks only the first n classes")
cert = json.load(open(sys.argv[1]))
a = cert['args']
lim = int(sys.argv[2]) if len(sys.argv) > 2 else a['classes']     # optional prefix check
print("regenerating: pool=%d orbit=%d walk=%d nu=%d sample=%d classes=%d scan=%d block=%d "
      "seed=%d" % (a['pool'], a['orbit'], a['walk'], a['nu'], a['sample'], lim, a['scan'],
                   a['block'], a['seed']), flush=True)
g7.select_backend(g7.device())
P = g7.spread_pool(a['pool'], walk=a['walk'], seed=a['seed'], orbit=a['orbit'],
                   nu=a['nu'], sample=a['sample'],
                   maxrows=a.get('maxrows', 24_000_000), verbose=False)
assert len(P) == cert['pool'], "pool size %d != certified %d" % (len(P), cert['pool'])
C, sizes, digests, worst, pdig = g7.build_classes(
    P, lim, seed=a['seed'], scan=a['scan'], block=a['block'], pchunk=a['pchunk'],
    runlen=a['runlen'], scanbig=a.get('scanbig') or None, nbig=a.get('nbig', 0),
    permute=a.get('permute', False), cachefree=a.get('cachefree', 0), verbose=False)
if cert.get('pool_digest') and cert['pool_digest'] != pdig:
    print("POOL DIGEST DIFFERS: certificate %s, regenerated %s" % (cert['pool_digest'][:16],
                                                                   pdig[:16]))
    print("The pool is assembled by torch primitives whose tie-breaking is unspecified, so it")
    print("is reproducible within one software stack and not guaranteed across them.  The")
    print("class sizes below will differ; that is the environment, not an error in the run.")
elif cert.get('pool_digest'):
    print("pool digest matches the certificate: %s" % pdig[:32])
else:
    print("the certificate predates pool digests; a mismatch below cannot be attributed")
h = hashlib.sha256()
for d in digests:
    h.update(d)
sz = list(sizes)
ok_sz = sz == cert['sizes'][:lim]
ok_h = (lim == a['classes']) and h.hexdigest() == cert['sha256']
print("pool %d, classes %d, sizes %d..%d, %d lines; worst |ip| = %d (<= 2), all distinct"
      % (len(P), len(C), min(sz), max(sz), sum(sz), worst))
print("sizes match: %s ; sha256 match: %s%s"
      % (ok_sz, ok_h, "" if lim == a['classes'] else "  (prefix of %d only)" % lim))
assert worst <= 2 and ok_sz and (ok_h or lim < a['classes'])
print("CERTIFICATE REPRODUCED")
