"""CPU ensemble of topology-changing 12D spherical packing experiments.

No unsuccessful numerical search is an upper bound or a rigidity certificate.
Each process uses one BLAS thread. Checkpoints include seeds and objective values.
"""
import os
for _key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS']:
    os.environ[_key] = '1'

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'


def normalize(x):
    return x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-100)


def metrics(x):
    x = normalize(x)
    g = x @ x.T
    values = g[np.triu_indices(len(x), 1)]
    excess = np.maximum(values - .5, 0)
    return dict(n=len(x), max_cos=float(values.max()),
                energy=float(excess @ excess),
                violations=int(np.count_nonzero(excess > 1e-12)),
                contacts_1e5=int(np.count_nonzero(values > .5-1e-5)))


class Hinge:
    def __init__(self, shape, target=.5, fixed=None, anchor=None, strength=0., weights=None):
        self.shape = shape
        self.target = target
        self.fixed = fixed
        self.anchor = anchor
        self.strength = strength
        self.weights = weights
        self.calls = 0

    def __call__(self, raw):
        self.calls += 1
        raw = raw.reshape(self.shape)
        norms = np.maximum(np.linalg.norm(raw, axis=1, keepdims=True), 1e-100)
        x = raw / norms
        g = np.maximum(x @ x.T - self.target, 0)
        np.fill_diagonal(g, 0)
        if self.weights is not None:
            value = .5 * np.sum(self.weights * g*g)
            grad = 2 * (self.weights*g) @ x
        else:
            value = .5 * np.sum(g*g)
            grad = 2 * g @ x
        if self.fixed is not None:
            cross = np.maximum(x @ self.fixed.T - self.target, 0)
            value += np.sum(cross*cross)
            grad += 2 * cross @ self.fixed
        if self.strength:
            dx = x - self.anchor
            value += self.strength * np.sum(dx*dx)
            grad += 2*self.strength*dx
        grad -= np.sum(grad*x, axis=1, keepdims=True)*x
        return float(value), (grad/norms).ravel()


def relax(x, steps, target=.5, fixed=None, anchor=None, strength=0., weights=None):
    fun = Hinge(x.shape, target, fixed, anchor, strength, weights)
    out = minimize(fun, normalize(x).ravel(), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=steps, maxfun=steps*3, maxcor=12,
                                ftol=1e-16, gtol=1e-10, maxls=30))
    return normalize(out.x.reshape(x.shape)), dict(calls=fun.calls, nit=out.nit,
                                                  objective=float(out.fun),
                                                  message=str(out.message))


def seed_configuration(route, rng, trial):
    canon = np.load(DATA/'canonical840.npy')
    record = np.load(DATA/'record841.npy')
    if route in ['quaternion26', 'record26', 'core25_control']:
        count = 25 if route == 'core25_control' else 26
        if route == 'quaternion26':
            core = canon[:816].copy()
            exceptional = np.zeros((count, 12))
            exceptional[:24] = canon[816:]
            exceptional[24:, 8:] = normalize(rng.normal(size=(count-24, 4)))
            exceptional[:, :8] += rng.normal(scale=.03, size=(count,8))
        else:
            core = record[:816].copy()
            exceptional = record[816:].copy()
            if count == 26:
                # Splitting different points explores different exceptional contact graphs.
                i = trial % 25
                z = exceptional[i].copy()
                tangent = rng.normal(size=12)
                tangent -= (tangent @ z)*z
                tangent /= np.linalg.norm(tangent)
                angle = [.25, .35, .45, .55][(trial//25) % 4]
                exceptional[i] = np.cos(angle)*z + np.sin(angle)*tangent
                exceptional = np.vstack([exceptional, np.cos(angle)*z-np.sin(angle)*tangent])
            exceptional += rng.normal(scale=[.0,.02,.06,.12][trial%4],size=exceptional.shape)
        return np.vstack([core, normalize(exceptional)]), 816
    if route == 'prune844':
        x = np.load(DATA/'published844.npy').copy()
        gram = x @ x.T
        np.fill_diagonal(gram, -1)
        bad = np.maximum(gram-.5,0)**2
        for _ in range(2):
            score = np.sum(bad,axis=1)
            if trial%3 == 0:
                i = int(rng.choice(len(x),p=score/score.sum()))
            elif trial%3 == 1:
                top = np.argsort(score)[-60:]
                i = int(rng.choice(top))
            else:
                i = int(rng.integers(len(x)))
            x = np.delete(x,i,axis=0)
            bad = np.delete(np.delete(bad,i,axis=0),i,axis=1)
        x += rng.normal(scale=[0.,.005,.02,.05][trial%4],size=x.shape)
        return normalize(x), 0
    if route == 'published842_control':
        return np.load(DATA/'published842.npy'), 0
    if route == 'record_global':
        x = record.copy()
        # Near the weak quaternionic factor, but without frozen core constraints.
        z = np.zeros(12)
        z[8:] = rng.normal(size=4)
        x = np.vstack([x, z/np.linalg.norm(z)])
        x += rng.normal(scale=[0.,.01,.03,.08][trial%4],size=x.shape)
        return normalize(x), 0
    raise ValueError(route)


def run_trial(spec):
    threadpool_limits(1)
    route, trial, seed, steps, rounds, outdir = spec
    outdir = Path(outdir)
    rng = np.random.default_rng(seed)
    start = time.monotonic()
    x, ncore = seed_configuration(route, rng, trial)
    initial = metrics(x)
    trace = []
    if ncore:
        # Allow the exceptional points to settle before letting the backbone respond.
        # This is only an initialization stage; the full core is then released.
        free, opt = relax(x[ncore:], min(steps,500), target=.5, fixed=x[:ncore])
        x[ncore:] = free
        trace.append(dict(stage='free_region', **metrics(x), **opt))
    best = None
    best_metric = None
    for r in range(rounds):
        x, opt = relax(x, steps)
        m = metrics(x)
        trace.append(dict(stage='full', round=r, **m, **opt))
        if best_metric is None or m['energy'] < best_metric['energy']:
            best = x.copy()
            best_metric = m
            np.save(outdir/f'{route}_{trial:04d}.npy',best)
        if m['max_cos'] < .5-1e-10:
            break
        if r+1 == rounds:
            break
        # Local collective perturbations can remove edges from a jammed contact graph.
        x = best.copy()
        if ncore and r % 3 != 2:
            count = int(rng.choice([4,8,16,26]))
            chosen = rng.choice(np.arange(ncore,len(x)),size=min(count,len(x)-ncore),replace=False)
            scale = float(rng.choice([.02,.05,.10,.20]))
            x[chosen] += rng.normal(scale=scale,size=(len(chosen),12))
            # Also release the most strongly interacting core points.
            score = np.maximum(x[:ncore] @ x[chosen].T - .46, 0).sum(axis=1)
            patch = np.argsort(score)[-int(rng.choice([16,32,64])):]
            x[patch] += rng.normal(scale=scale*.3,size=(len(patch),12))
        else:
            scale = float(rng.choice([.003,.01,.03]))
            x += rng.normal(scale=scale,size=x.shape)
        x = normalize(x)
    result = dict(route=route,trial=trial,seed=seed,steps=steps,rounds=rounds,
                  initial=initial,best=best_metric,trace=trace,
                  elapsed_seconds=time.monotonic()-start)
    (outdir/f'{route}_{trial:04d}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--routes',nargs='+',default=['quaternion26','record26','prune844','record_global'])
    parser.add_argument('--trials',type=int,default=16)
    parser.add_argument('--workers',type=int,default=32)
    parser.add_argument('--steps',type=int,default=1500)
    parser.add_argument('--rounds',type=int,default=3)
    parser.add_argument('--seed',type=int,default=20261104)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    specs = [(route,t,args.seed+100000*j+t,args.steps,args.rounds,str(args.out.resolve()))
             for t in range(args.trials) for j,route in enumerate(args.routes)]
    config = vars(args).copy()
    config['out'] = str(args.out.resolve())
    config['started_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    config['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run_trial,s):s for s in specs}
        for future in as_completed(futures):
            result=future.result()
            results.append(result)
            print(json.dumps(dict(route=result['route'],trial=result['trial'],
                                  best=result['best'],seconds=round(result['elapsed_seconds'],2))),flush=True)
            (args.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__ == '__main__':
    main()
