#!/usr/bin/env python3
"""The 32000-class family is REGENERABLE: check that, on a prefix, end to end.

    python verify_rebuild.py [N]        # default 60 classes, about a minute

Dimensions 81-95 read the sizes of a family of 32 000 pairwise-disjoint Gamma_72 classes --
47 150 230 lines, several gigabytes -- so the family itself is not shipped; only
`data/class_sizes_32000.npy` is, and the README says the classes are rebuilt deterministically
from `data/class_best_2106.npy` and `data/gamma72_gens.npy`.

That claim went untested, and it was false: `build_classes.py` and `materialise.py` were copied
from the working tree with its directory layout and looked for their inputs in a `dim73-80/`
directory that this package does not have, so neither would run at all.  Nothing exercised
them, so nothing noticed.  This script exercises them:

  1. build_classes.py rebuilds the first N classes as (automorphism, survivor mask);
  2. materialise.py expands that to explicit coefficient vectors;
  3. verify_classes.py's checks are applied to every one of them -- norm 8, every off-diagonal
     |<u,v>| <= 2, all lines distinct;
  4. the rebuilt sizes are compared with the first N entries of the shipped size table.

Step 4 is the one that ties the shipped numbers to the shipped scripts.  It uses a scratch
directory and leaves the package unchanged.

HOW DEEP IS DEEP ENOUGH.  N is a prefix in BUILD order, but final81-95.py reads the T+P
LARGEST classes and data/class_sizes_32000.npy is not sorted, so a prefix of T+P does NOT
cover the dimension that reads T+P.  The depth each dimension really needs is the largest
build index among its T+P largest classes, and it is printed at the end of this run.  It is
about 1.2 to 1.4 times T+P: dimension 85 reads 302 classes but needs 447 rebuilt, and only
N = 32 000 covers dimension 95.
"""
import io
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60

print(__doc__)
tmp = tempfile.mkdtemp(prefix='kn_rebuild_')
env = dict(os.environ, KN_STATE=os.path.join(tmp, 'state'), KN_OUT=os.path.join(tmp, 'classes'))
try:
    for step, args in (('build_classes.py', [str(N)]), ('materialise.py', [])):
        r = subprocess.run([sys.executable, os.path.join(HERE, step)] + args,
                           env=env, capture_output=True, text=True, cwd=HERE)
        if r.returncode != 0:
            print(r.stdout[-2000:]); print(r.stderr[-2000:])
            raise SystemExit("%s failed -- the family is NOT regenerable from what is shipped"
                             % step)
        print("   %-20s ok" % step)

    batches = sorted(f for f in os.listdir(env['KN_OUT']) if f.endswith('.npz'))
    assert batches, "materialise.py produced nothing"
    sys.path.insert(0, HERE)
    import verify_classes as vc

    G = np.load(os.path.join(PKG, 'data', 'gamma72_gram.npy')).astype(np.int64)
    total, keys, sizes = 0, [], []
    rng = np.random.default_rng(20260823)
    r = rng.integers(-(2 ** 51), 2 ** 51, size=72, dtype=np.int64)
    for b in batches:
        for name, C in vc.load(os.path.join(env['KN_OUT'], b)):
            A = C.astype(np.int64)
            assert A.shape[1] == 72, (name, A.shape)
            nrm = np.einsum('ij,jk,ik->i', A, G, A)
            assert (nrm == 8).all(), "%s: a vector is not minimal" % name
            IP = A @ G @ A.T
            np.fill_diagonal(IP, 0)
            assert int(np.abs(IP).max()) <= 2, "%s: not a 60-degree code" % name
            keys.append(np.abs(A @ r))
            sizes.append(len(A))
            total += len(A)
    K = np.concatenate(keys)
    assert len(np.unique(K)) == total, "the rebuilt classes are not pairwise disjoint"
    print("   %d classes, %d lines: every norm 8, every |off-diagonal| <= 2, all distinct"
          % (len(sizes), total))

    ship = np.load(os.path.join(PKG, 'data', 'class_sizes_32000.npy')).ravel().astype(np.int64)
    got = np.array(sizes, dtype=np.int64)
    assert len(ship) == 32000, "the shipped size table is not 32000 long"
    assert (got == ship[:len(got)]).all(), (
        "the rebuilt sizes differ from the shipped table:\n  rebuilt %s\n  shipped %s"
        % (got[:12].tolist(), ship[:len(got)][:12].tolist()))
    print("   rebuilt sizes agree with data/class_sizes_32000.npy on all %d" % len(got))
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
print("REGENERABLE: the shipped scripts reproduce the shipped size table, and every class they")
print("produce is a genuine 60-degree code of Gamma_72 minimal lines, disjoint from the rest.")

# ---------------------------------------------------------------------------------------
# What this N covered, and what each dimension would need.  Derived from the shipped size
# table and rows81-95.json, so it cannot drift from either.
import json

_rows = os.path.join(PKG, 'rows81-95.json')
_sz = os.path.join(PKG, 'data', 'class_sizes_32000.npy')
if os.path.exists(_rows) and os.path.exists(_sz):
    _a = np.load(_sz).ravel().astype(np.int64)
    _ord = np.argsort(-_a, kind='stable')
    print()
    print("   what a prefix of N covers, on the image family (the GPU family is separate)")
    print("   dim   classes read   deepest build index among them   covered by N = %d" % N)
    for _r in json.load(io.open(_rows, encoding='utf-8')):
        _n = _r['T'] + _r['P']
        if _n > len(_a):
            print("   %3d %10d   -- more classes than this family has" % (_r['dim'], _n))
            continue
        _d = int(_ord[:_n].max()) + 1
        print("   %3d %10d %25d   %s" % (_r['dim'], _n, _d, 'yes' if N >= _d else 'no'))
    print()
    print("   A dimension is covered only when N reaches its depth: the driver reads the")
    print("   LARGEST classes, this script rebuilds the FIRST ones, and they are not the")
    print("   same set.  Dimensions 86-95 take the GPU family in the shipped result, which")
    print("   this script does not rebuild; on the image family alone they are still")
    print("   records, by the margin README.md states.")
