"""Derive the registered totals from pinned parent masks and signed deltas."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
PARENT=HERE.parent/'dim49-60-p48-pointwise-refill'
POLES=(2,6,12,24,40,72,126,240,306,510,604,841)
CLASSES=(1,3,6,12,20,36,63,120,153,255,302,420)

def results():
 metadata=json.loads((HERE/'certificate.json').read_text())
 if [r['dimension'] for r in metadata['rows']]!=list(range(49,61)):raise ValueError('Dimension registry mismatch')
 out=[]
 for row,poles,classes in zip(metadata['rows'],POLES,CLASSES):
  d=row['dimension'];parent=PARENT/f'heads{d}.npz';delta=HERE/f'p48_signed_delta{d}.json.gz'
  if hashlib.sha256(parent.read_bytes()).hexdigest()!=row['parent_sha256']:raise ValueError('Pinned parent mismatch')
  if hashlib.sha256(delta.read_bytes()).hexdigest()!=row['delta_sha256']:raise ValueError('Signed delta mismatch')
  with np.load(parent,allow_pickle=False) as a:
   bits=a['retained_bits'];extra=a['extra_points'];offsets=a['extra_offsets']
   if bits.dtype!=np.uint8 or bits.shape!=(classes,973):raise ValueError('Parent mask shape')
   mask=np.unpackbits(bits,axis=1)
   if mask[:,7783:].any() or extra.ndim!=2 or extra.shape[1]!=48:raise ValueError('Parent mask or extra coordinates')
   if offsets.shape!=(classes+1,) or offsets[0]!=0 or offsets[-1]!=len(extra) or np.any(np.diff(offsets)<0):raise ValueError('Parent extra partition')
   old_signed=2*(int(mask[:,:7783].sum())+len(extra))
  decoded=gzip.decompress(delta.read_bytes())
  if hashlib.sha256(decoded).hexdigest()!=row['decoded_delta_sha256']:raise ValueError('Decoded signed delta mismatch')
  C=json.loads(decoded);removed=sum(len(c['remove']) for c in C['changes']);added=sum(len(c['add']) for c in C['changes'])
  total=52416000+old_signed-removed+added+poles
  if total!=row['lower_bound'] or total!=C['bound'] or added-removed!=row['gain']:raise ValueError('Derived total disagreement')
  if total<=row['previous_repository_bound']:raise ValueError('No improvement')
  out.append(dict(dimension=d,total=total))
 return out

if __name__=='__main__':
 for row in results():print(json.dumps(row))
