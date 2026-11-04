"""Validate saved evidence, exercise SAT controls, aggregate and hash the run."""
import hashlib,json,os
from types import SimpleNamespace
import numpy as np
from pysat.solvers import Solver
from model import HERE,DATA,OLD,N,labels_from_holes,verify_labels
from certify_trade_barrier import parameters,tuples
from certify_halves import verify
from solve_sat import build

all_holes=np.concatenate([np.load(DATA/'recovery.npz')['holes'],
                          np.load(DATA/'nonlinear_quotient_holes.npy')])
for H in all_holes:
    parameters(H)
    assert verify_labels(labels_from_holes(H))==960
assert not tuples(12)
records=json.loads((DATA/'half_certificates.json').read_text())['certificates']
verify(records)
saved=[]
for path in sorted(DATA.glob('*_labels.npy')):
    L=np.load(path);n=verify_labels(L)
    assert n<=960, 'Potential improvement must receive full geometric certification'
    saved.append(dict(file=path.name,selected=n,status='PASS'))

# Full-label CNF: a valid supplied configuration must be SAT; corrupting one
# label to point to a selected word must be UNSAT. Auxiliaries remain free.
args=SimpleNamespace(nonlinear=False,initial=0,mode='labels',cap=64,radius=None,
                     fixed_half=None,half_cuts=False)
cnf,H,x=build(args);L=labels_from_holes(H)
def assumptions(labels):
    return ([u+1 if H[u] else -u-1 for u in range(1024)] +
            [int(x[u,p])*(1 if labels[u]==p else -1)
             for u in range(1024) for p in range(21)])
with Solver(name='cadical195',bootstrap_with=cnf) as solver:
    assert solver.solve(assumptions=assumptions(L))
    corrupt=L.copy()
    for u in np.flatnonzero(L>=0):
        directions=np.flatnonzero(L[N[u]]>=0)
        if len(directions):corrupt[u]=directions[0];break
    assert not solver.solve(assumptions=assumptions(corrupt))

pending=[];runs=[]
for log in sorted(DATA.glob('*.log')):
    path=DATA/(log.stem+'_result.json')
    if not path.exists():pending.append(log.stem)
    else:
        r=json.loads(path.read_text())
        if 'run' in r:runs.append(r)
assert not pending,pending
recovery=json.loads((DATA/'recovery.json').read_text())
geometry=json.loads((DATA/'geometry_verification.json').read_text())
assert geometry['status']=='PASS' and geometry['count']==31560
summary=dict(dimension=21,verified_lower_bound=31560,improvement_found=False,
             recovery=dict(successes=32,tests=32,maximum_seconds=recovery['maximum_seconds'],
                           scope=recovery['assumed_structure']),
             nonlinear_seeds=8,geometry=geometry,
             trade_barrier=json.loads((DATA/'trade_barrier_certificate.json').read_text()),
             rational_halves=json.loads((DATA/'halves_verification.json').read_text()),
             full_sat_positive_control='PASS',full_sat_invalid_label_control='PASS',
             saved_labels=saved,pending_runs=pending,
             status_counts={status:sum(r['status']==status for r in runs)
                            for status in sorted({r['status'] for r in runs})},
             runs=runs)
(DATA/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
paths=list(HERE.glob('*.py'))+[HERE/'README.md',HERE/'.gitignore']
paths += [p for p in DATA.iterdir() if p.suffix in ('.json','.npz','.npy','.log')
          and p.name!='manifest.json']
paths += [OLD/'component21.npz',OLD/'baseline21.npz']
manifest={os.path.relpath(p,HERE).replace(os.sep,'/'):
          hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
(DATA/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(status='PASS',verified_lower_bound=31560,
                      recovery_successes=32,nonlinear_seeds=8,
                      status_counts=summary['status_counts'],pending_runs=pending,
                      hashed_files=len(manifest))),flush=True)
