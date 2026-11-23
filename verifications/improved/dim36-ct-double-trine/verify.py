#!/usr/bin/env python3
import argparse,hashlib,importlib.util,itertools,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--previous',type=Path,default=HERE.parent/'dim36-ct-triangle-refill')
 args=p.parse_args()
 src=args.previous.resolve()
 W=json.loads((HERE/'witness.json').read_text())
 for f,h in W['sha256'].items():
  d=(src/f).read_bytes();assert hashlib.sha256(d).hexdigest()==h,(f,'unexpected source')
 lib=importlib.util.spec_from_file_location('source36',src/'verify.py')
 m=importlib.util.module_from_spec(lib);sys.modules['source36']=m;lib.loader.exec_module(m)
 g=m.read_geometry();assert g.points==502902
 labels=np.flatnonzero(g.minima==4)
 R=g.roots[labels].reshape(-1,12)
 assert np.array_equal((R@g.gram*R).sum(1),np.full(len(R),4))
 enc={int(a):i for i,a in enumerate(labels)}
 old=sorted(g.final)
 I=np.array([[3*enc[a]+b for a,b in zip(row[:3],row[3:])]for row in old],dtype=np.int64)
 assert len(old)==99198
 inserted=set()
 for c in W['insert_full_triples']:
  assert len(c)==3 and all(int(y) in enc for y in c)
  for t in itertools.product(range(3),repeat=3):
   row=tuple(c)+t;assert row not in inserted
   inserted.add(row)
 assert not inserted&g.final
 J=np.array([[3*enc[a]+b for a,b in zip(row[:3],row[3:])]for row in sorted(inserted)],dtype=np.int64)
 assert len(J)==27*len(W['insert_full_triples'])
 T=R@g.gram@R.T
 newdot=sum(T[J[:,i,None],J[:,i]]for i in range(3))
 assert np.all(np.diag(newdot)==12)
 assert np.max(newdot-12*np.eye(len(J),dtype=np.int64))<=6
 # Shape is checked by ordinary integer dot, not a discovery conflict graph.
 D=set()
 mx=-13
 for j0 in range(0,len(J),24):
  q=J[j0:j0+24]
  d=sum(T[q[:,i,None],I[:,i]]for i in range(3))
  mask=(d>6).any(axis=0)
  D.update(old[i] for i in np.flatnonzero(mask))
  if not mask.all():mx=max(mx,int(d[:,~mask].max()))
 expected=set(map(tuple,W['remove_rows']))
 assert D==expected and len(D)==len(W['remove_rows']) and D<=g.final
 # The predecessor already audited existing rows; check every new-unchanged base pair.
 crosses=m.verify_cross(g,R,J)
 result=g.points-len(D)+len(J)
 assert result==W['points'] and result>g.points
 print(json.dumps(dict(status='PASS',dimension=36,previous=g.points,points=result,added=len(J),removed=len(D),
                        exact_new_old_root_products=len(J)*len(I),new_base_family_checks=crosses),indent=2))

if __name__=='__main__':
 assert __debug__
 main()
