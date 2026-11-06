"""Derive registered populations from the retained ownership masks."""
from pathlib import Path
import json
import numpy as np
D=Path(__file__).resolve().parent

def results():
 rows=[]
 for r in json.loads((D/'results.json').read_text()):
  a=np.load(D/f'classes_{r["k"]}.npz',allow_pickle=False);b=np.unpackbits(a['retained_bits'],axis=1)
  assert len(b)==r['k'] and not b[:,7783:].any() and np.array_equal(b[:,:7783].sum(1),a['sizes'])
  n=52416000+2*int(a['sizes'].sum())+2*r['k'];assert n==r['new_bound']
  rows.append({'dimension':r['dimension'],'total':n})
 return rows
if __name__=='__main__':
 for row in results():print(json.dumps(row))
