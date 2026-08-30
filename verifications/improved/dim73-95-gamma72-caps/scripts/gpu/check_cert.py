#!/usr/bin/env python3
"""Check the GPU class-family certificate without a GPU.

    python check_cert.py                 (from this directory)

The family of pairwise-disjoint Gamma_72 classes used by dimensions 81-95 is certified by a
SEED, not by its 6.5 GB of lines: every random decision in the pipeline is made by numpy's
PCG64 and the device does nothing but exact integer GEMM, so `regenerate.py cert_prod.json`
rebuilds the identical family anywhere and compares the SHA-256.  That needs a GPU to be
quick; this script checks everything that does not:

  * the sizes array actually used by final81-95.py is the one in the certificate;
  * the measured |<u,v>| histogram of the pool is consistent with the EXACT distribution of
    the Gamma_72 shell, which follows from its minimal vectors being a spherical 11-design;
  * the reported total and the top-T sums are what the sizes give.
"""
import io
import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'data'))

SHELL_LINES = {0: 1301829375, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
SHELL_TOTAL = sum(SHELL_LINES.values())
SHELL_P = (SHELL_LINES[3] + SHELL_LINES[4]) / SHELL_TOTAL

cert = json.load(open(os.path.join(HERE, 'cert_prod.json')))
sz = np.load(os.path.join(DATA, 'class_sizes_gpu.npy'))
print(__doc__)
assert list(sz) == cert['sizes'], "class_sizes_gpu.npy does not match the certificate"
print("sizes match the certificate: %d classes, %d..%d, %d lines in total"
      % (len(sz), sz.min(), sz.max(), sz.sum()))
assert int(sz.sum()) == cert['total_lines']

h, tot = cert.get('ip_histogram'), cert.get('ip_pairs')
if h:
    print()
    print("  |ip|      measured        exact shell      ratio")
    for c in range(5):
        print("   %d     %.7f      %.7f      %.4f"
              % (c, h[str(c)] / tot, SHELL_LINES[c] / SHELL_TOTAL,
                 (h[str(c)] / tot) / (SHELL_LINES[c] / SHELL_TOTAL)))
    p = sum(h[str(c)] for c in range(3, 9)) / float(tot)
    assert abs(p - cert['p_conflict']) < 1e-9
else:
    # A certificate written from a resumable state (cert_from_state.py) carries the measured
    # density but not the raw histogram, which only run.py computes.
    p = cert['p_conflict']
print("conflict density %.6f against the shell's exact %.7f (%+.1f%%)"
      % (p, SHELL_P, 100 * (p / SHELL_P - 1)))
assert p >= SHELL_P * 0.98, "measured density below the exact shell value -- impossible"

if cert.get('partial'):
    print("PARTIAL run: %d of the %d classes asked for; dimensions needing more than %d "
          "classes keep the other family" % (cert['classes_built'], cert['classes_requested'],
                                             cert['classes_built']))
    assert cert['alive_lines'] + cert['total_lines'] == cert['pool'], \
        "disjointness invariant violated"
    print("   invariant: %d alive + %d used = %d pool lines"
          % (cert['alive_lines'], cert['total_lines'], cert['pool']))

srt = np.sort(sz)[::-1].astype(np.int64)
print()
print("  greedy predicts ln(1 + scan*p)/p = %.0f lines a class; realised mean %.0f"
      % (np.log(1 + cert['args']['scan'] * p) / p, sz.mean()))
# What dimension 95 actually reads is T triples and P pairs, worth 4 and 2 lines each --
# not "4x the top 31050", which was a literal left behind by an earlier partition and was
# both the wrong count of classes (31 072) and the wrong weight.  Read from the driver's own
# output so it cannot drift again.
_rows = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'rows81-95.json')
if os.path.exists(_rows):
    _r95 = [r for r in json.load(io.open(_rows, encoding='utf-8')) if r['dim'] == 95][0]
    _T, _P = _r95['T'], _r95['P']
    print("  dimension 95 reads %d classes as %d triples and %d pairs: 4*%d + 2*%d + %d "
          "poles = %d" % (_T + _P, _T, _P, int(srt[:_T].sum()),
                          int(srt[_T:_T + _P].sum()), _r95['poles'],
                          4 * int(srt[:_T].sum()) + 2 * int(srt[_T:_T + _P].sum())
                          + _r95['poles']))
else:
    print("  rows81-95.json is not here, so what dimension 95 reads cannot be shown; "
          "run final81-95.py first")
print("  seed %d, pool %d, scan %d/%d, block %d, backend %s"
      % (cert['args']['seed'], cert['pool'], cert['args']['scan'],
         cert['args'].get('scanbig', 0), cert['args']['block'], cert['backend']))
print("  sha256 of the classes: %s" % cert['sha256'])
print()
print("ALL CHECKS PASSED (run regenerate.py on a GPU to rebuild and match the sha256)")
