"""Three candidate-specific corruptions must fail both head checkers."""
from pathlib import Path
import json,struct,subprocess,tempfile
import numpy as np
from verify_single_class import check,HERE
z=np.load(HERE/'p48_signed_class_lns.npz');X=z['points'];G=z['gram'];keys=set(map(tuple,X))
# This preserves all norms and distinctness while reversing a real negative-3 pair.
bad=None
for i,x in enumerate(X):
 if tuple(-x) in keys:continue
 h=X@G@x
 if np.any(h==-3):
  bad=X.copy();bad[i]=-x;assert len(np.unique(bad,axis=0))==len(X);break
assert bad is not None
cases={'duplicate':np.r_[X,X[:1]],'wrong_norm':X.copy(),'positive_product_three':bad}
cases['wrong_norm'][0,0]+=1
assert cases['wrong_norm'][0]@G@cases['wrong_norm'][0]!=6
results={}
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp);exe=tmp/'checker'
 subprocess.run(['g++','-O3','-std=c++17',str(HERE/'check_signed_class.cpp'),'-o',str(exe)],check=True)
 for name,a in cases.items():
  p=tmp/(name+'.npz');np.savez_compressed(p,points=a,gram=G)
  rejected=False
  try:check(p)
  except (AssertionError,ValueError):rejected=True
  assert rejected
  binary=tmp/(name+'.bin');binary.write_bytes(struct.pack('<Q',len(a))+G.astype('<i8').tobytes()+a.astype('<i8').tobytes())
  result=subprocess.run([str(exe),str(binary)],capture_output=True,text=True)
  assert result.returncode!=0
  results[name]={'python':'REJECTED','cpp':'REJECTED','cpp_exit':result.returncode,'reason':result.stderr.strip()}
print(json.dumps({'status':'PASS','controls':results},indent=2))
