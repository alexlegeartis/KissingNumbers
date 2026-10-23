#!/usr/bin/env python3
"""Search half-lattice cosets of Gamma_72; candidates and certificates are separate.

Coordinates are rows in the catalogue basis, with Gram matrix G.  A parity p
denotes p+2 Gamma_72.  Representatives of norm r give heads of norm r/4.
The enumerator uses a rational LDL decomposition, enclosed by outward rounded
IEEE double intervals.  Its input also retains the exact integral Gram matrix.

Examples:
    python gamma72_cosets.py prepare
    python gamma72_cosets.py enumerate --seed 1 --radius 24 --seconds 30
"""
import argparse
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'dim73-80'
PREFIX = HERE / 'gamma72_cosets'


def gram_lll(gram, delta=.99):
    """LLL of an exact Gram matrix; floating GSO only selects integral moves.

    The returned change of basis and Gram identity are checked exactly.  LLL
    quality has no bearing on correctness of the later exhaustive enumeration.
    """
    g = gram.astype(np.int64).copy()
    b = np.eye(len(g), dtype=np.int64)
    n = len(g)
    chol = np.linalg.cholesky(g.astype(float))
    d = np.diag(chol)**2
    mu = chol / np.diag(chol)[None, :]
    k, swaps, reductions = 1, 0, 0
    while k < n:
        for j in range(k-1, -1, -1):
            q = int(round(mu[k, j]))
            if q:
                b[k] -= q*b[j]
                old = g[k].copy()
                g[k] -= q*g[j]
                g[:, k] -= q*g[:, j]
                g[k, k] = old[k] - 2*q*old[j] + q*q*g[j, j]
                mu[k, :j+1] -= q*mu[j, :j+1]
                reductions += 1
        if d[k] >= (delta-mu[k, k-1]**2)*d[k-1]:
            k += 1
        else:
            g[[k, k-1]] = g[[k-1, k]]
            g[:, [k, k-1]] = g[:, [k-1, k]]
            b[[k, k-1]] = b[[k-1, k]]
            chol = np.linalg.cholesky(g.astype(float))
            d = np.diag(chol)**2
            mu = chol / np.diag(chol)[None, :]
            swaps += 1
            k = max(k-1, 1)
    assert np.array_equal(b @ gram @ b.T, g)
    return g, b, {'swaps': swaps, 'reductions': reductions,
                  'gs_min': float(d.min()), 'gs_max': float(d.max())}


def rational_ldl(g):
    n = len(g)
    mu = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    d = []
    for i in range(n):
        d.append(Fraction(int(g[i, i])) - sum(mu[i][k]**2*d[k] for k in range(i)))
        assert d[-1] > 0
        for j in range(i+1, n):
            mu[j][i] = (int(g[j, i])-sum(mu[j][k]*mu[i][k]*d[k]
                                        for k in range(i))) / d[i]
    return d, mu


def enclose(q):
    f = float(q)
    lo = f if Fraction.from_float(f) <= q else math.nextafter(f, -math.inf)
    hi = f if Fraction.from_float(f) >= q else math.nextafter(f, math.inf)
    assert Fraction.from_float(lo) <= q <= Fraction.from_float(hi)
    return lo, hi


def write_ldl(g, data):
    d, mu = rational_ldl(g)
    with data.open('w') as f:
        f.write(str(len(g))+'\n')
        for row in g:
            f.write(' '.join(map(str, row))+'\n')
        for x in d:
            f.write('%.17g %.17g\n' % enclose(x))
        for row in mu:
            for x in row:
                f.write('%.17g %.17g\n' % enclose(x))


def prepare():
    original = np.load(SOURCE / 'gamma72_gram.npy').astype(np.int64)
    start = time.time()
    g, b, stats = gram_lll(original)
    write_ldl(g, HERE / 'gamma72_cosets_ldl.txt')
    np.savez(HERE / 'gamma72_cosets_basis.npz', gram=g, basis=b)
    stats['seconds'] = time.time()-start
    print(json.dumps(stats), flush=True)
    return stats


def enumerate_coset(args):
    if not (HERE / 'gamma72_cosets_basis.npz').exists():
        prepare()
    dat = np.load(HERE / 'gamma72_cosets_basis.npz')
    rng = np.random.default_rng(args.seed)
    parity = rng.integers(0, 2, size=72)
    if args.parity:
        parity = np.array([int(c) for c in args.parity], dtype=np.int64)
        assert len(parity) == 72
    inp = HERE / ('gamma72_coset_%d.input.txt' % args.seed)
    inp.write_text(' '.join(map(str, parity))+'\n')
    out = HERE / ('gamma72_coset_%d' % args.seed)
    cmd = [str(HERE / 'gamma72_cosets_enum.exe'),
           str(HERE / 'gamma72_cosets_ldl.txt'), str(inp),
           str(args.radius), str(args.seconds), str(out)]
    env = dict(os.environ)
    env['PATH'] = 'C:/msys64/mingw64/bin;'+env.get('PATH', '')
    subprocess.run(cmd, check=True, env=env)
    if not Path(str(out)+'.vectors.txt').stat().st_size:
        print('No representatives found before completion/time limit.')
        return
    rows = np.loadtxt(str(out)+'.vectors.txt', dtype=np.int64, ndmin=2)
    if rows.size:
        assert np.all((rows-parity) % 2 == 0)
        norms = np.einsum('ij,jk,ik->i', rows, dat['gram'], rows)
        assert norms.max() <= args.radius
        unique, counts = np.unique(norms, return_counts=True)
        print('exact checked representatives:', dict(zip(map(int, unique), map(int, counts))))
        np.savez(str(out)+'.npz', parity=parity, reduced_reps=rows,
                 catalogue_reps=rows @ dat['basis'], norms=norms)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['prepare', 'enumerate'])
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--radius', type=int, default=24)
    ap.add_argument('--seconds', type=float, default=30)
    ap.add_argument('--parity')
    args = ap.parse_args()
    if args.command == 'prepare':
        prepare()
    else:
        enumerate_coset(args)


if __name__ == '__main__':
    main()
