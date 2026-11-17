from pathlib import Path
import numpy as np,sympy as sp,json,hashlib,subprocess
from fractions import Fraction as Fr
P=Path(__file__).resolve().parent.parent;A=P/'independent_audit';P0=P.parent/'mixed_shell71_20261117';import sys
dim=int(sys.argv[1]); embedding={68:16,69:25}[dim]; expected={68:'5ca894e106eb934eebe7c8605fae99846a3cfcec22f43954ba8a06027f1d53f9',69:'d05673d3241427cbfcbc80c77dbf3e6dab33fd466170c64029061901a3ebe3d2'}[dim]
raw=(P/f'scan_best{dim}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==expected;target=json.loads(raw)
X=np.load(P/f'kernel_{embedding}.npz');F=X['F'];N=X['normals'];K=X['K'];H=X['H'];rows=X['solve_rows'].tolist();tail=[i for i in range(24)if i not in rows];M=F@N.T
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
 z=subprocess.run([str(A/'count128'),str(path),str(P0/'fiber8_table.bin'),'unused','8'],capture_output=True,text=True,check=True)
 return json.loads(z.stdout)
report=[]
for i,z in enumerate(target['fibers']):
 h=z['h'];t=Fr(z.get('height','0'));r=8-Fr(sum(x*x for x in h),16);aa=(8-t*t)/r;assert aa>0
 for w in target['fibers'][:i+1]:
  k=w['h'];u=Fr(w.get('height','0'));rr=8-Fr(sum(x*x for x in k),16);bb=(8-u*u)/rr;dot=4-Fr(sum(x*y for x,y in zip(h,k)),16);c=4-t*u;assert c>=0 and c*c>=aa*bb*dot*dot
 j=R@h;assert np.all(j%4==0);j=j//4;z8=count(j,8);pure=count(j//2,4)['hist'].get('4',0)if np.all(j%2==0)else 0;value=z8['weighted']+pure;assert value>=z['count'];report.append(dict(h=h,claimed=z['count'],replayed=value,pure=pure,enumeration=z8));print(i,value,flush=True)
assert target.get('poles',0)==0 and sum(z['count']for z in target['fibers'])==target['bound']
assert len({tuple(z['h']) for z in target['fibers']})==len(target['fibers'])
result=dict(status='PASS',sha256=hashlib.sha256(raw).hexdigest(),dimension=dim,bound=target['bound'],pair_checks=len(report)*(len(report)+1)//2,rows=report,scope='Exact rational geometry and saturated affine mapping; direct repeated traversal with exact integer leaves provides lower counts, not completeness. Inherits audited Leech completion weights and published Gamma72 minimum8.')
(A/f'report_scan{dim}.json').write_text(json.dumps(result,indent=2));print('PASS',result['bound'],flush=True)
