#!/usr/bin/env python3
"""Bounded CP-SAT search on an exact weighted orbit graph; every incumbent is saved."""
import argparse
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from ortools.sat.python import cp_model
from orbits17 import DATA, read_graph


def solve(name, seconds, workers, seed, require_gain):
    started = time.monotonic()
    z = np.load(DATA / (name+'.npz'))
    weights = z['weights'].astype(int)
    divisor = math.gcd(*map(int, weights))
    w = weights // divisor
    off, nb = read_graph(DATA / (name+'.graph'))
    n = len(weights)
    model = cp_model.CpModel()
    x = [model.new_bool_var(f'x{i}') for i in range(n)]
    cover = DATA / (name+'.cliques')
    if cover.exists():
        for line in cover.open():
            clique = list(map(int,line.split()))
            assert clique[0] == len(clique)-1
            model.add_at_most_one(x[i] for i in clique[1:])
    elif int(z['dimension']) == 20 and 'dense_representatives' in z:
        from prepare20 import basis_of, reduction
        import itertools
        _, coordinates, _ = reduction(basis_of(z['subgroup']), 20)
        ball = [sum(1 << p for p in positions) for k in range(3)
                for positions in itertools.combinations(range(20), k)]
        ball = coordinates(ball)
        assert len(np.unique(ball)) == 211
        ndense = len(z['dense_representatives'])
        # Every distance<=4 difference splits into two weight<=2 errors.
        covered = set(map(int, np.bitwise_xor.outer(ball,ball).ravel()))-{0}
        actual = set(map(int, nb[off[0]:off[1]][nb[off[0]:off[1]] < ndense]))
        assert actual == covered
        for i in range(ndense):
            model.add_at_most_one(x[int(j)] for j in i ^ ball)
            sparse = nb[off[i]:off[i+1]]
            sparse = sparse[sparse >= ndense]
            model.add_bool_and(~x[int(j)] for j in sparse).only_enforce_if(x[i])
        for i in range(ndense,n):
            for j in nb[off[i]:off[i+1]]:
                if j > i:
                    model.add_at_most_one(x[i],x[int(j)])
    else:
        for i in range(n):
            for j in nb[off[i]:off[i+1]]:
                if j > i:
                    model.add_at_most_one(x[i], x[int(j)])
    objective = sum(int(wi)*xi for wi, xi in zip(w, x))
    initial = set(map(int, z['seed']))
    baseline = sum(int(w[i]) for i in initial)
    fixed_vertex = None
    if int(z['dimension']) == 20:
        # Every all-sign translation is an automorphism of the full 20D pool.
        # A solution beating the root+octad base must contain a dense point.
        fixed_vertex = 0
    elif 'representatives' in z:
        # Deliberately search the branch containing a flat point. No unrestricted
        # upper bound is inferred for a hypothetical solution without any flats.
        reps = z['representatives']
        fixed_vertex = next(i for i in initial if 61440 <= int(reps[i]) < 63488)
    if fixed_vertex is not None:
        model.add(x[fixed_vertex] == 1)
    model.add(objective >= baseline + int(require_gain))
    model.maximize(objective)
    for i in range(n):
        model.add_hint(x[i], int(i in initial))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.log_search_progress = False
    solver.parameters.cp_model_probing_level = 0
    solver.parameters.cp_model_presolve = False
    solver.parameters.symmetry_level = 0
    solver.parameters.linearization_level = 0
    fixed = int(z['fixed'])
    run = f'{name}_s{seed}' + ('_gain' if require_gain else '')

    class Save(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.best = -1
        def on_solution_callback(self):
            score = int(round(self.objective_value))
            if score <= self.best:
                return
            self.best = score
            selected = np.asarray([i for i in range(n) if self.value(x[i])], dtype=np.int32)
            # Check against the independent, saved graph, not the solver model.
            chosen = np.zeros(n, dtype=bool)
            chosen[selected] = True
            assert all(not chosen[nb[off[i]:off[i+1]]].any() for i in selected)
            np.save(DATA / (run+'_best.npy'), selected)
            print(json.dumps(dict(event='incumbent', run=run, count=fixed+divisor*score,
                                  seconds=round(time.monotonic()-started,3))), flush=True)

    print(json.dumps(dict(event='solve_start', run=run, variables=n,
                          edges=int(off[-1])//2, model_seconds=time.monotonic()-started,
                          workers=workers, limit=seconds, baseline=fixed+divisor*baseline)), flush=True)
    cb = Save()
    status = solver.solve(model, cb)
    has_solution = status in (cp_model.OPTIMAL,cp_model.FEASIBLE)
    bound = fixed+divisor*int(math.floor(solver.best_objective_bound+1e-6))
    report = dict(run=run, status=solver.status_name(status), dimension=int(z['dimension']),
                  best_count=(fixed+divisor*int(round(solver.objective_value)) if has_solution else fixed+divisor*baseline),
                  restricted_upper_bound=(bound if bound >= fixed+divisor*baseline else None),
                  baseline=fixed+divisor*baseline, solver_seconds=solver.wall_time,
                  total_seconds=time.monotonic()-started, workers=workers,
                  conflicts=solver.num_conflicts, branches=solver.num_branches,
                  fixed_vertex=fixed_vertex,
                  scope='Only this subgroup-invariant finite candidate pool; not an unrestricted kissing upper bound.')
    (DATA / (run+'_result.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('name')
    p.add_argument('--seconds', type=float, default=600)
    p.add_argument('--workers', type=int, default=16)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--require-gain', action='store_true')
    a = p.parse_args()
    solve(a.name,a.seconds,a.workers,a.seed,a.require_gain)
