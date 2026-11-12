"""Isolated MILP worker, so the parent can enforce a real wall deadline."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz
from scipy.optimize import Bounds, LinearConstraint, milp

ap = argparse.ArgumentParser()
ap.add_argument('--directory', type=Path, required=True)
ap.add_argument('--seconds', type=float, default=90)
args = ap.parse_args()
f = np.load(args.directory/'model_data.npz', allow_pickle=False)
c, ub = f['c'], f['ub']
A = load_npz(args.directory/'model_constraints.npz')
result = milp(c, integrality=np.ones(len(c)), bounds=Bounds(np.zeros(len(c)), np.ones(len(c))),
              constraints=LinearConstraint(A, -np.inf, ub),
              options={'time_limit': args.seconds, 'threads': 2})
report = dict(solver_status=int(result.status), solver_message=str(result.message),
              has_incumbent=result.x is not None,
              objective=None if result.fun is None else float(result.fun))
if result.x is not None:
    np.save(args.directory/'solver_solution.npy', result.x)
(args.directory/'solver_result.json').write_text(json.dumps(report, indent=2)+'\n')
