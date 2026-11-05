"""Reproduce a 13D mixed-norm exchange search; no record is assumed."""
import argparse, json, pathlib, sys
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from discovery.core import RunContext
from discovery.engines.projective_exchange import search, certify

p = argparse.ArgumentParser()
p.add_argument('--seconds', type=float, default=120)
p.add_argument('--gain', type=int, default=33)
p.add_argument('--output', type=pathlib.Path)
a = p.parse_args()
z = np.load(pathlib.Path(__file__).with_name('pool.npz'))
X, Q, normal = (z[k].astype(np.int64) for k in ('source','candidates','normal'))
assert X.shape == (1122,14) and Q.shape == (6224,14)
assert np.any(normal) and not np.any(X @ normal) and not np.any(Q @ normal)
ctx = RunContext(ROOT, {'max_dimension': 96, 'max_materialized_points': 20000,
    'max_codewords': 10000, 'max_pairs': 30000000, 'max_rank_work': 500000000,
    'max_gram_entries': 100000000, 'solver_seconds': max(a.seconds,1)})
proposal = search([], {'source':X.tolist(), 'candidates':Q.tolist(),
                       'minimum_gain':a.gain, 'seconds':a.seconds}, ctx)
code = certify([proposal], {}, ctx)
print(json.dumps({'source_points':len(X), 'candidates':len(Q),
  'minimum_gain':a.gain, 'solver_status':proposal.data['solver_status'],
  'selected':len(proposal.data['selected_indices']),
  'removed':len(proposal.data['removed_indices']), 'exact_facts':code.facts}, indent=2))
if a.output and proposal.data['selected_indices']:
    a.output.write_text(json.dumps(proposal.data))
