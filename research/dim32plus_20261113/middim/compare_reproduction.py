"""Compare regenerated discovery outputs with the frozen exact construction.
Only elapsed-time metadata is ignored. No geometric computation is repeated.
"""
from pathlib import Path
import gzip,json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CANONICAL=ROOT/'verifications/improved/dim49-60-p48-signed-refill'

def main():
 if not __debug__:raise SystemExit('Assertions must be enabled.')
 with np.load(HERE/'p48_signed_class_lns.npz',allow_pickle=False) as a, np.load(CANONICAL/'p48_signed_class_lns.npz',allow_pickle=False) as b:
  assert set(a.files)==set(b.files)=={'points','gram'}
  assert np.array_equal(a['gram'],b['gram'])
  assert len(a['points'])==len(b['points'])
  assert set(map(tuple,a['points']))==set(map(tuple,b['points']))
  points=len(a['points'])
 dimensions=[]
 for d in range(49,61):
  name=f'p48_signed_delta{d}.json.gz'
  local=json.loads(gzip.decompress((HERE/name).read_bytes()))
  frozen=json.loads(gzip.decompress((CANONICAL/name).read_bytes()))
  local.pop('seconds',None);frozen.pop('seconds',None)
  assert local==frozen,('Regenerated signed delta differs',d)
  dimensions.append(d)
 print(json.dumps(dict(status='PASS',signed_class_points=points,identical_constructions=dimensions,ignored_metadata=['elapsed seconds only']),indent=2))

if __name__=='__main__':main()
