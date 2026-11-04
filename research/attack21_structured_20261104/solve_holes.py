"""Independent domination or necessary domination relaxation, with local branching.

No orbit invariance is imposed. A radius bounds the number of baseline holes
removed; when absent all 1024 hole indicators are free except translation.
"""
import argparse, json, time
import numpy as np
from ortools.sat.python import cp_model
from model import C, N, SIDE, DATA, labels_from_holes, verify_component_geometry
from certify_trade_barrier import parameters, tuples


def solve(a):
    started = time.monotonic()
    initial = (np.load(DATA / 'nonlinear_quotient_holes.npy')[a.initial] if a.nonlinear
               else np.load(DATA / 'recovery.npz')['holes'][a.initial])
    m = cp_model.CpModel()
    h = [m.new_bool_var(f'h{u}') for u in range(1024)]
    gain = m.new_bool_var('gain')
    m.add(sum(h)<=63).only_enforce_if(gain)
    m.add(sum(h)>=64).only_enforce_if(~gain)
    for u in range(1024):
        m.add_bool_or([h[u]] + [h[v] for v in N[u]])
        if not a.relax:
            for v in N[u]:
                if v > u: m.add_at_most_one(h[u], h[v])
    m.add(sum(h) <= a.cap)
    if a.fixed_half is not None or a.half_cuts:
        from certify_halves import quotient, verify
        if a.fixed_half is not None:
            _,T=quotient(a.fixed_half)
            for u in np.searchsorted(C,T): m.add(h[u]==1)
        if a.half_cuts:
            records=json.loads((DATA/'half_certificates.json').read_text())['certificates']
            verify(records)
            kernel=np.load(DATA/'recovery.npz')['holes'][0]
            for record in records:
                inside=np.zeros(1024,dtype=bool)
                inside[np.searchsorted(C,record['fixed_subgroup'])]=True
                for half in (inside,kernel&~inside):
                    m.add_bool_or([~h[u] for u in np.flatnonzero(half)]).only_enforce_if(gain)
    # Bipartite covering inequalities are redundant but expose useful arithmetic.
    sides = [sum(h[u] for u in np.flatnonzero(SIDE == s)) for s in (0, 1)]
    m.add(sides[0] + 21*sides[1] >= 512)
    m.add(sides[1] + 21*sides[0] >= 512)
    if a.radius is not None:
        m.add(sum(~h[u] for u in np.flatnonzero(initial)) <= a.radius)
    else:
        m.add(h[0] == 1)  # translate any nonempty hole set to contain zero
    if a.trade_cuts:
        parameters(initial)
        counts = []
        for old in (True, False):
            for side in (0,1):
                ids = np.flatnonzero((initial == old) & (SIDE == side))
                v = m.new_int_var(0,len(ids),f'count_{old}_{side}')
                m.add(v == sum((~h[u] if old else h[u]) for u in ids))
                counts.append(v)
        ra,rb,sa,sb = counts
        m.add(12*ra <= 17*sa+sb); m.add(12*rb <= 17*sb+sa)
        m.add(ra <= 3*sb); m.add(rb <= 3*sa)
        m.add(ra+rb >= 13).only_enforce_if(gain)
        if a.radius is not None:
            m.add_allowed_assignments(counts,tuples(a.radius)).only_enforce_if(gain)
    if a.hint:
        for u in range(1024): m.add_hint(h[u], int(initial[u]))
    if a.optimize: m.minimize(sum(h))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = a.seconds
    s.parameters.num_search_workers = a.workers
    s.parameters.random_seed = a.seed
    s.parameters.log_search_progress = True
    s.parameters.linearization_level = a.linearization
    if a.no_presolve: s.parameters.cp_model_presolve = False
    s.parameters.symmetry_level = a.symmetry
    class Save(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            H = np.array([self.value(x) for x in h], dtype=bool)
            np.save(DATA / (a.name+'_holes.npy'), H)
            if not a.relax:
                L = labels_from_holes(H)
                np.save(DATA / (a.name+'_labels.npy'), L)
            print(json.dumps(dict(event='solution', run=a.name, holes=int(H.sum()),
                                  seconds=time.monotonic()-started)), flush=True)
    status = s.solve(m, Save())
    report = dict(run=a.name, status=s.status_name(status), cap=a.cap,
                  independent=not a.relax, radius=a.radius, initial=a.initial, nonlinear=a.nonlinear,
                  hint=a.hint, optimize=a.optimize, workers=a.workers, trade_cuts=a.trade_cuts,
                  fixed_half=a.fixed_half,half_cuts=a.half_cuts,
                  seconds=s.wall_time, total_seconds=time.monotonic()-started,
                  branches=s.num_branches, conflicts=s.num_conflicts,
                  response_stats=s.response_stats())
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        H = np.array([s.value(x) for x in h], dtype=bool)
        report['holes'] = int(H.sum())
        if not a.relax:
            report['geometry'] = verify_component_geometry(labels_from_holes(H))
    if a.optimize: report['hole_lower_bound'] = s.best_objective_bound
    (DATA / (a.name+'_result.json')).write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--name', required=True)
    p.add_argument('--cap', type=int, default=63)
    p.add_argument('--seconds', type=float, default=300)
    p.add_argument('--workers', type=int, default=16)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--initial', type=int, default=0)
    p.add_argument('--nonlinear', action='store_true')
    p.add_argument('--radius', type=int)
    p.add_argument('--relax', action='store_true')
    p.add_argument('--hint', action='store_true')
    p.add_argument('--trade-cuts', action='store_true')
    p.add_argument('--fixed-half',type=int)
    p.add_argument('--half-cuts',action='store_true')
    p.add_argument('--optimize', action='store_true')
    p.add_argument('--no-presolve', action='store_true')
    p.add_argument('--linearization', type=int, default=2)
    p.add_argument('--symmetry', type=int, default=2)
    solve(p.parse_args())
