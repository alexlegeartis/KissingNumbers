"""Portable, read-only exact replay of six P48 coverage certificates."""
from pathlib import Path
import hashlib,itertools,json
import numpy as np
D=Path(__file__).resolve().parent

def sign(a,b,r):
 a,b,r=int(a),int(b),int(r)
 if not b or not r:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if (a>0)==(b>0):return 1 if a>0 else -1
 q=a*a-r*b*b
 return ((q>0)-(q<0))*(1 if a>0 else -1)

def tails():
 roots=[]
 for i,j in itertools.combinations(range(8),2):
  for s,t in itertools.product((-2,2),repeat=2):
   x=np.zeros(8,dtype=np.int64);x[i]=s;x[j]=t;roots.append(x)
 for x in itertools.product((-1,1),repeat=8):
  if sum(v<0 for v in x)%2==0:roots.append(x)
 roots=np.array(roots,dtype=np.int64)
 for k,n in ((6,72),(7,126),(8,240)):
  x=roots
  if k<8:x=x[x[:,6]+x[:,7]==0]
  if k<7:x=x[x[:,5]+x[:,7]==0]
  assert len(x)==n and len(set(map(tuple,x)))==n
  assert all(tuple(-v) in set(map(tuple,x)) for v in x)
  h=x@x.T;assert np.all(np.diag(h)==8);np.fill_diagonal(h,0);assert h.max()<=4
  # The displayed independent linear equations put E6/E7 in at most k dimensions.
 for k,n in ((9,306),(10,510),(11,604)):
  z=np.load(D/'source'/f'capdirs_{k}.npz',allow_pickle=False)
  A=z['A'].astype(np.int64);B=z['B'].astype(np.int64);c=z['c'].astype(np.int64);r=int(z['R'].item());N=int(z['norm'].item())
  assert A.shape==B.shape==(n,k) and c.shape==(k,) and np.all(c>0) and N>0
  bound=k*int(c.max())*(int(abs(A).max())+int(abs(B).max()))**2*max(1,r)
  assert 4*bound<2**63
  U=(A*c)@A.T+r*((B*c)@B.T);V=(A*c)@B.T+(B*c)@A.T
  assert np.all(np.diag(U)==N) and not np.diag(V).any()
  keys={tuple(np.r_[a,b]) for a,b in zip(A,B)};assert len(keys)==n
  assert all(tuple(-np.r_[a,b]) in keys for a,b in zip(A,B))
  for i in range(n):
   for j in range(i):assert sign(2*U[i,j]-N,2*V[i,j],r)<=0

def main():
 if not __debug__:raise SystemExit('Run without -O: certificate assertions must be enabled.')
 m=json.loads((D/'manifest.json').read_text())
 for name,h in m['files'].items():assert hashlib.sha256((D/name).read_bytes()).hexdigest()==h,name
 tails()
 z=np.load(D/'source/family_refilled.npz',allow_pickle=False);base=z['base'].astype(np.int64);G=z['gram'].astype(np.int64)
 assert base.shape==(7783,48) and G.shape==(48,48) and np.array_equal(G,G.T)
 assert abs(base).max()<=17 and abs(G).max()<=10
 for start in range(0,len(base),256):
  h=base[start:start+256]@G@base.T
  for j,row in enumerate(h):
   i=start+j;assert row[i]==6;assert np.all(abs(row[:i])<=2) and np.all(abs(row[i+1:])<=2)
 out=[]
 for r in json.loads((D/'results.json').read_text()):
  a=np.load(D/f'classes_{r["k"]}.npz',allow_pickle=False);T=a['transforms'].astype(np.int64);bits=np.unpackbits(a['retained_bits'],axis=1)
  assert T.shape==(r['k'],48,48) and abs(T).max()<=21
  assert not bits[:,7783:].any();bits=bits[:,:7783].astype(bool)
  assert np.array_equal(bits.sum(1),a['sizes']) and int(a['sizes'].sum())==r['final']
  assert int(z['sizes'][:r['k']].sum())==r['initial']
  # Coarse bounds: TGT^t <=48^2*21^2*10; x=base*T <=48*17*21.
  # All products and partial sums, including xGx^t, are safely below 2^63.
  assert 48**2*10*(48*17*21)**2<2**63
  allrows=[];kept=[]
  for t,mask in zip(T,bits):
   assert np.array_equal(t@G@t.T,G)
   x=base@t;j=np.argmax(x!=0,axis=1);x*=np.sign(x[np.arange(len(x)),j])[:,None]
   assert np.max(abs(x))<128
   x=x.astype(np.int8);allrows.append(x);kept.append(x[mask])
  X=np.concatenate(allrows);keys=X.view(np.dtype((np.void,48))).ravel();u,first=np.unique(keys,return_index=True)
  expected=np.zeros(len(X),bool);expected[first]=True;assert np.array_equal(expected.reshape(bits.shape),bits)
  K=np.concatenate(kept).view(np.dtype((np.void,48))).ravel();assert len(np.unique(K))==len(K)==len(u)==r['final']
  poles=2*r['k'];N=52416000+2*len(K)+poles;assert N==r['new_bound']
  assert N-(52416000+2*r['initial']+poles)==r['point_gain']
  out.append(dict(dimension=r['dimension'],total=N,gain=r['point_gain'],selected_lines=len(K),status='PASS'))
 print(json.dumps({'status':'PASS','rows':out,'base_pairs':30283653,'scope':'Exact finite heads, directions and poles; published P48 minimum six and 52416000 shell population remain inherited.'},indent=2))
if __name__=='__main__':main()
