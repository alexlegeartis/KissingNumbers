#!/usr/bin/env python3
"""Optimize nonlinear line labels on one 1024-word component, with optional subgroups."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from prepare21 import DATA,eta
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'compute17_20_20261104'))
from prepare20 import basis_of


def solve(rank,seconds,workers,seed,gain):
    started=time.monotonic();z=np.load(DATA/'component21.npz')
    component=z['component'];fullN=z['neighbors'];labels=z['labels']
    baseline_words=np.load(DATA/'baseline21.npz')['words']
    kernel=[int(m) for m in component if eta(int(baseline_words[m]))==0]
    rng=np.random.default_rng(seed);basis=basis_of(rng.permutation(kernel))[:rank]
    H=[0]
    for b in basis:H += [m^b for m in H]
    assert len(H)==2**rank
    full_index={int(m):i for i,m in enumerate(component)}
    orbit=np.full(1024,-1,dtype=np.int32);reps=[]
    for u,m in enumerate(component):
        if orbit[u]>=0:continue
        members=[full_index[int(m)^h] for h in H]
        assert np.all(orbit[members]<0)
        orbit[members]=len(reps);reps.append(u)
        assert len(set(map(int,labels[members])))==1
    N=orbit[fullN[reps]];n=len(reps);initial=labels[reps]
    assert np.all(N!=np.arange(n)[:,None])
    model=cp_model.CpModel();hole=[model.new_bool_var(f'z{u}') for u in range(n)]
    color=[[model.new_bool_var(f'x{u}_{p}') for p in range(21)] for u in range(n)]
    for u in range(n):
        model.add(sum(color[u])+hole[u]==1)
        for p in range(21):
            model.add_implication(color[u][p],hole[int(N[u,p])])
            neighbors=set(map(int,N[u]));neighbors.discard(int(N[u,p]))
            model.add_bool_and(~color[v][p] for v in neighbors).only_enforce_if(color[u][p])
        # Explicit domination cuts make the relevant small-cover obstruction visible.
        model.add(hole[u]+sum(hole[v] for v in set(map(int,N[u])))>=1)
        model.add_hint(hole[u],int(initial[u]<0))
        for p in range(21):model.add_hint(color[u][p],int(initial[u]==p))
    baseline=int(np.sum(initial<0))
    model.add(sum(hole)<=baseline-int(gain))
    model.minimize(sum(hole))
    # Every admissible selection has a hole. Translation preserves the line rule,
    # so one hole may be placed at zero without loss, also in the orbit models.
    model.add(hole[int(orbit[full_index[0]])]==1)
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=workers;solver.parameters.random_seed=seed
    solver.parameters.symmetry_level=0;solver.parameters.cp_model_probing_level=0
    solver.parameters.cp_model_presolve=False;solver.parameters.linearization_level=0
    run=f'cp21_r{rank}_s{seed}'+('_gain' if gain else '')
    class Save(cp_model.CpSolverSolutionCallback):
        def __init__(self):super().__init__();self.best=baseline+1
        def on_solution_callback(self):
            h=int(round(self.objective_value))
            if h>=self.best:return
            self.best=h
            assigned=np.asarray([next((p for p in range(21) if self.value(color[u][p])),-1) for u in range(n)],dtype=np.int32)[orbit]
            for u,p in enumerate(assigned):
                if p<0:continue
                assert assigned[fullN[u,p]]<0
                assert all(assigned[v]!=p for v in fullN[u] if assigned[v]>=0)
            selected=int(np.sum(assigned>=0));total=27720+4*selected
            np.save(DATA/(run+'_labels.npy'),assigned)
            print(json.dumps(dict(event='incumbent',run=run,component_selected=selected,total_if_replicated=total,seconds=time.monotonic()-started)),flush=True)
            if selected>960:self.stop_search()
    print(json.dumps(dict(event='start',run=run,words=n,rank=rank,variables=22*n,workers=workers,seconds=seconds)),flush=True)
    cb=Save();status=solver.solve(model,cb)
    report=dict(run=run,status=solver.status_name(status),rank=rank,workers=workers,
                seconds=solver.wall_time,total_seconds=time.monotonic()-started,
                best_component_size=1024-int(round(solver.objective_value))*2**rank if status in (cp_model.OPTIMAL,cp_model.FEASIBLE) else None,
                hole_lower_bound=int(np.ceil(solver.best_objective_bound))*2**rank,
                branches=solver.num_branches,conflicts=solver.num_conflicts,
                target='more than 960 points in one component; replicate in all four for a 21D improvement')
    (DATA/(run+'_result.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rank',type=int,default=0);p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--workers',type=int,default=16);p.add_argument('--seed',type=int,default=1);p.add_argument('--gain',action='store_true')
    a=p.parse_args();solve(a.rank,a.seconds,a.workers,a.seed,a.gain)
