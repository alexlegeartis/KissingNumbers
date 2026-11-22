"""Exact finer block evaluation of the pinned Gamma72 independence recurrence."""
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
import gzip,hashlib,importlib.util,json,os,sys
from math import isqrt,ceil
from pathlib import Path
HERE=Path(__file__).resolve().parent
P=Path(os.getenv('ROOTED_PARENT',str(HERE.parent/'dim74-85-gamma72-rooted-exchange'))).resolve()
PINS={'verify.py':'652a0c1f0295a9c58a922b8f087d2000390e5b00a2d14c6dc420f22a66bf6cce','witness.json':'b57e1284104a3130f6ef63bca34de4014ed39cd507e36a715f22066ec2b8ba6a','results.json':'8590cdd78a51c12e8df7b3335a97897e47354def2aadc18042784465c8c32125','oriented.json.gz':'e238a0566472d40bfcc54f47f2008fe57cb3197692419d9feebd3f4a3f79f32d','oriented.py':'5f4ed8f637659e192181b7ca9e68f0ef2cc2d7a2d4a370f31c4d9a7eb1348f54'}
scale=10**32
N=15281788354560;N8=6218175600
@lru_cache(None)
def lf(D,T,m):
 a0=isqrt(2*T)+2
 while a0*(a0-1)<=2*T:a0+=1
 assert a0<=D
 W=max(1,(D-a0+m)//m)
 g=0
 for a in range(a0,D+1,W):
  b=min(D,a+W-1);L=b-a+1
  den=a*a+1-2*T;n=den-a-1;cd=a*den;cn=2*(a+T)
  assert n>=0 and 0<cn<cd
  gc=max(0,cd-L*cn);bc=scale*L*n*cd;bd=b*den*(cd+L*cn)
  g=(gc*g*bd+bc*cd)//(cd*bd)
 return F(scale+g,scale)
def replay():
 for f,h in PINS.items():assert hashlib.sha256((P/f).read_bytes()).hexdigest()==h,f
 os.environ['HARMONIC_PARENT']=os.getenv('HARMONIC_PARENT',str(P.parent/'dim74-85-gamma72-harmonic-latitudes'))
 sys.path.insert(0,str(P))
 sp=importlib.util.spec_from_file_location("rooted_parent_verifier",P/'verify.py')
 q=importlib.util.module_from_spec(sp);sp.loader.exec_module(q)
 old=q.replay()
 assert old==json.loads((P/'results.json').read_text())
 dims={d['n']:d for d in json.loads((P/'witness.json').read_text())['results']}
 with gzip.open(P/'oriented.json.gz','rt')as f:items={int(z['d']):z for z in json.load(f)}
 R=[]
 for r in old:
  d=r['dimension'];k=d-72;T=r['triangles'];base=r['selected']
  if d>=77:
   M,D=int(dims[d]['M']),int(dims[d]['D'])
   assert 2*ceil(F(N//2*M,D+1)*lf(D,T,512))==base
   L=2*ceil(F(N//2*M,D+1)*lf(D,T,32768))
  else:
   D=Counter(items[d]['degrees'])
   assert ceil(sum(F(N*n,v+1)*lf(v,T,512)for v,n in D.items()))==base
   L=ceil(sum(F(N*n,v+1)*lf(v,T,32768)for v,n in D.items()))
  gain=L-base
  assert gain>0
  R.append({'dimension':d,'previous':r['lower_bound'],'lower_bound':r['lower_bound']+gain,'gain':gain,'selected':L,'triangles':T})
  print(d,R[-1]['lower_bound'],gain,flush=True)
 return R
if __name__=='__main__':
 assert __debug__
 ans=replay()
 if '--record' in sys.argv:(HERE/'results.json').write_text(json.dumps(ans,indent=2)+'\n')
 else:assert ans==json.loads((HERE/'results.json').read_text())
