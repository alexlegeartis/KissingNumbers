"""Search a fixed exterior for off-alphabet sphere insertions in 26D.

The reduced coordinates are a chosen Leech octad plus the two tail axes.
Every exterior point is represented by its actual projection; deduplication
preserves all fixed-exterior inequalities. Numerical outputs are candidates,
never lower-bound certificates.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


def normalize(x):
    return x / np.linalg.norm(x, axis=1)[:, None]


def objective(flat, m, d, boundary):
    raw = flat.reshape(m, d)
    norms = np.linalg.norm(raw, axis=1)
    p = raw / norms[:, None]
    v = np.maximum(p @ p.T - 0.5, 0)
    np.fill_diagonal(v, 0)
    w = np.maximum(p @ boundary.T - 0.5, 0)
    value = 0.5 * np.sum(v * v) + np.sum(w * w)
    gradient = 2 * (v @ p + w @ boundary)
    gradient -= np.sum(gradient * p, axis=1)[:, None] * p
    return value, (gradient / norms[:, None]).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, default=Path(__file__).with_name('normal26.npz'))
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--trials', type=int, default=64)
    ap.add_argument('--free-all', action='store_true')
    ap.add_argument('--seed', type=int, default=20261113)
    ap.add_argument('--maxiter', type=int, default=1500)
    a = ap.parse_args()
    z = np.load(a.input)
    y = z['Y'] / np.sqrt(32)
    b = z['B'] / np.sqrt(32)
    if a.free_all:
        moving = y
    else:
        moving = y[364:]
        b = np.vstack([b, y[:364]])
    rng = np.random.default_rng(a.seed)
    rows = []
    best = float('inf')
    start = time.monotonic()
    for trial in range(a.trials):
        seed = moving.copy()
        i = int(rng.integers(364 if a.free_all else 0, len(moving)))
        v = rng.normal(size=seed.shape[1])
        v -= (v @ seed[i]) * seed[i]
        v /= np.linalg.norm(v)
        theta = rng.uniform(0.3, 0.6)
        new = np.cos(theta)*seed[i] + np.sin(theta)*v
        seed[i] = np.cos(theta)*seed[i] - np.sin(theta)*v
        seed = np.vstack([seed, new])
        if trial % 4 == 3:
            seed = normalize(seed + rng.normal(size=seed.shape)*0.04)
        m, d = seed.shape
        res = minimize(objective, seed.ravel(), args=(m, d, b), jac=True,
                       method='L-BFGS-B', options={'maxiter': a.maxiter,
                       'ftol': 1e-16, 'gtol': 2e-11, 'maxls':40,'maxcor':30})
        p = normalize(res.x.reshape(m, d))
        g = p @ p.T
        np.fill_diagonal(g, -1)
        violation = max(g.max()-0.5, (p@b.T).max()-0.5)
        row = dict(trial=trial, loss=float(res.fun), violation=float(violation),
                   iterations=int(res.nit), seconds=time.monotonic()-start)
        rows.append(row)
        if res.fun < best:
            best = float(res.fun)
            np.savez_compressed(a.out.with_suffix('.npz'), candidate=p*np.sqrt(32),
                                fixed_boundary=b*np.sqrt(32), indices=z['indices'])
            print(json.dumps(dict(row, best=True)), flush=True)
        a.out.write_text(json.dumps(dict(scope='Numerical finite-restart search; no certificate',
                            seed=a.seed, free_all=a.free_all, moving=len(moving),
                            candidate_count=m,boundary_count=len(b),trials=rows), indent=2)+'\n')
        if violation < 1e-10:
            print('NUMERICAL FEASIBILITY REQUIRES EXACTIFICATION', flush=True)
            break


if __name__ == '__main__':
    main()
