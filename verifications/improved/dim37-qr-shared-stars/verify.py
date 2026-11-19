"""Exact 37D star insertion on the 533,864-point QR incumbent."""
from itertools import combinations,product
from pathlib import Path
import numpy as np,json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"dim37-38-qr-merge"))
from verify37 import replay as base_replay
from supports37 import supports
P=Path(__file__).resolve().parent
def stars():
 d=json.loads((P/'star_recipe.json').read_text())
 assert d['special']==36
 H=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]],dtype='i8')
 H=np.r_[H,-H]
 def fibre(c,leaves):
  assert len(set(leaves))==4 and c not in leaves and 0<=c<37 and min(leaves)>=0 and max(leaves)<37
  for s in (-4,4):
   X=np.zeros((8,37),'i8');X[:,c]=s;X[:,leaves]=2*H
   yield X
 out=[]
 for b in d['partition']:out.extend(fibre(36,b))
 for c,b in d['ordinary'].items():out.extend(fibre(int(c),[*b,36]))
 return np.concatenate(out)
def replay():
 original=base_replay()
 assert original['points']==533864
 X=stars();n=len(X);assert n==720 and np.all((X*X).sum(1)==32)
 assert len(np.unique(X,axis=0))==n
 assert max(np.sort(abs(X),axis=1)[:,-1])<5
 P0=X@X.T;np.fill_diagonal(P0,-32);assert P0.max()<=16
 S=supports();J=np.array([[j for j in range(36)if (s>>j)&1]for s in S],dtype='i8')
 # An absolute-sum bound covers all 128 signings in each incumbent octad.
 for i in range(0,n,64):
  assert (2*abs(X[i:i+64])[:,J].sum(2)).max()<=16
 # Incumbent dense coordinates are regular magnitude two and marked
 # magnitudes four or one. Use worst sign alignment in each norm class.
 for a in (1,4):
  D=2*abs(X[:,:36]).sum(1)+a*abs(X[:,36]);assert 4*D.max()**2<=32*(144+a*a)
 R=[]
 for i,j in combinations(range(37),2):
  for u,v in product((-4,4),repeat=2):
   r=np.zeros(37,'i8');r[i]=u;r[j]=v;R.append(r)
 R=np.array(R);assert len(R)==2664
 B=X@R.T>16
 cut=B.any(0); assert cut.sum()==576
 roots=R[~cut];assert (roots@X.T).max()<=16
 points=original['points']+n-int(cut.sum())
 assert points==534008
 return dict(status='PASS',dimension=37,points=points,baseline=original['points'],stars=n,deleted_roots=int(cut.sum()),retained_roots=len(roots),star_star_max=int(P0.max()))
if __name__=='__main__':
 import verify_six
 base=replay()
 extra=verify_six.replay_six(stars())
 print(json.dumps({**base,**extra,'points':extra['total']},sort_keys=True))
