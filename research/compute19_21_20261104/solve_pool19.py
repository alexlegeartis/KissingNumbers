#!/usr/bin/env python3
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from baselines import DATA
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'compute17_20_20261104'))
from orbits17 import read_graph

p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=180);p.add_argument('--workers',type=int,default=16);p.add_argument('--exact',action='store_true');a=p.parse_args()
name='polar19_exact' if a.exact else 'polar19'
started=time.monotonic();off,nb=read_graph(DATA/(name+'.graph'));n=len(off)-1
model=cp_model.CpModel();x=[model.new_bool_var(f'x{i}') for i in range(n)]
for i in range(n):
    for j in nb[off[i]:off[i+1]]:
        if j>i:model.add_at_most_one(x[i],x[int(j)])
model.add(sum(x)>=1793);model.maximize(sum(x))
for i in range(n):model.add_hint(x[i],int(i<1792))
solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.num_search_workers=a.workers
solver.parameters.symmetry_level=0;solver.parameters.cp_model_probing_level=0;solver.parameters.linearization_level=0
class Save(cp_model.CpSolverSolutionCallback):
    def on_solution_callback(self):
        selected=np.asarray([i for i in range(n) if self.value(x[i])],dtype=np.int32)
        assert len(selected)>1792
        chosen=set(map(int,selected));assert all(chosen.isdisjoint(map(int,nb[off[i]:off[i+1]])) for i in selected)
        np.save(DATA/('cp19_'+name+'_best.npy'),selected)
        print(json.dumps(dict(event='candidate',count=10668+len(selected),seconds=time.monotonic()-started)),flush=True)
        self.stop_search()
print(json.dumps(dict(event='start',vertices=n,edges=int(off[-1])//2,seconds=a.seconds,workers=a.workers)),flush=True)
status=solver.solve(model,Save())
report=dict(status=solver.status_name(status),seconds=solver.wall_time,total_seconds=time.monotonic()-started,
            best_count=10668+int(solver.objective_value) if status in (cp_model.OPTIMAL,cp_model.FEASIBLE) else 12460,
            branches=solver.num_branches,conflicts=solver.num_conflicts,
            exact_field_arithmetic=a.exact,
            scope='Sampled rational polar directions and the published auxiliary layer, with the complete base fixed.')
(DATA/('cp19_'+name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
