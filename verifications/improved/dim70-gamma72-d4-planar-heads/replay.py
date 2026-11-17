from pathlib import Path
import numpy as np,sympy as sp,json,hashlib,subprocess
from fractions import Fraction as Fr
import sys
P=Path(__file__).resolve().parent; A=Path(sys.argv[1]); P0=P.parent/'dim68-91-lattice-code-continuation/mixed_shell71_20261117'; target=json.load(open(P/'candidate.json'))
X=np.load(P/'kernel_0.npz');F=X['F'];N=X['normals'];K=X['K'];H=X['H'];rows=X['solve_rows'].tolist();tail=[i for i in range(24)if i not in rows];M=F@N.T
assert np.array_equal(F,np.load(P0/'leech_structure.npz')['gram'])
assert abs(sp.Matrix(M[rows]).det())==1 and abs(sp.Matrix(K[:,tail]).det())==1 and np.array_equal(K@M,np.zeros((20,4),int))
assert np.array_equal(H,K@F@K.T)
inv=np.array(sp.Matrix(M[rows]).inv(),np.int64);Hi=np.array(sp.Matrix(H).inv(),float);L,D=sp.Matrix(H).LDLdecomposition(hermitian=False);L=np.array(L,float);d=np.diag(np.array(D,float));keys=[sum((int(x)&1)<<j for j,x in enumerate(row))for row in K]
R=np.array([[0,2,-2,0],[2,-2,0,0],[0,0,2,-2],[0,0,2,2]],np.int64);assert np.array_equal(2*N@F@N.T,R@R.T)
def count(j,radius):
 v=np.zeros(24,np.int64);v[rows]=np.array(j)@inv;assert np.array_equal(v@M,j)
 h=K@F@v;q=int(v@F@v);sh=Hi@h;c=float(q-h@sh);key=sum((int(x)&1)<<k for k,x in enumerate(v));path=A/'affine.input'
 with open(path,'w')as f:
  f.write(f'20 {radius} {q} {c:.17g} {key}\n')
  for row in H:f.write(' '.join(map(str,row))+'\n')
  for arr in[h,keys,d,sh]:f.write(' '.join(map(str,arr))+'\n')
  for row in L:f.write(' '.join(map(str,row))+'\n')
 z=subprocess.run([str(A/'count128'),str(path),str(A/'fiber8_table.bin'),'unused','8'],capture_output=True,text=True,check=True)
 return json.loads(z.stdout)

assert len({tuple(z['h']) for z in target})==len(target)
report=[]
for i,z in enumerate(target):
 h=z['h']; u=z['u']; v=z['v']; r=128-sum(t*t for t in h); T=3*r+u*u+3*v*v
 assert r>0 and 0<T<=384
 for w in target[:i]:
  k=w['h']; S=3*(128-sum(t*t for t in k))+w['u']**2+3*w['v']**2
  B=3*(64-sum(x*y for x,y in zip(h,k)))+u*w['u']+3*v*w['v']
  assert B<=0 or 4*B*B<=T*S
 j=R@h;assert np.all(j%4==0);j=j//4
 z8=count(j,8); pure=count(j//2,4)['hist'].get('4',0) if np.all(j%2==0) else 0
 value=z8['weighted']+pure;assert value>=z['count'],(h,value,z['count'])
 report.append(dict(h=h,claimed=z['count'],replayed=value,pure=pure,enumeration=z8))
 (A/'partial.json').write_text(json.dumps(report,indent=2));print(i,value,flush=True)
assert sum(z['count'] for z in target)==1361134720
(A/'RESULT.json').write_text(json.dumps(dict(status='PASS',bound=1361134720,rows=report),indent=2))
print('PASS',flush=True)
