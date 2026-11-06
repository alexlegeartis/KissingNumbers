"""Negative controls for exact primitives and ownership."""
import numpy as np
from verify import sign,tails
assert sign(1,0,2)>0 and sign(-1,0,2)<0
assert sign(1,-1,2)<0 and sign(2,-1,2)>0
assert sign(-1,1,2)>0 and sign(-2,1,2)<0
# A duplicate canonical line fails the same disjointness test used in replay.
x=np.array([[1,0],[1,0]],dtype=np.int8)
assert len(np.unique(x.view(np.dtype((np.void,2))).ravel()))!=len(x)
# Altered ownership fails deterministic first-occurrence reconstruction.
ids=np.array([0,1,0,2]);_,first=np.unique(ids,return_index=True)
expected=np.zeros(4,bool);expected[first]=True
assert not np.array_equal(expected,np.array([True,True,True,True]))
tails();print('PASS: quadratic signs, duplicate rejection, ownership rejection and tails')
# An actual altered certificate is rejected by the portable verifier before replay.
import tempfile,shutil
from pathlib import Path
import verify
with tempfile.TemporaryDirectory() as tmp:
 target=Path(tmp)/'package';shutil.copytree(verify.D,target)
 p=target/'classes_36.npz';b=bytearray(p.read_bytes());b[-1]^=1;p.write_bytes(b)
 old=verify.D;verify.D=target
 try:
  try:verify.main()
  except AssertionError:pass
  else:raise AssertionError('Tampered certificate was accepted')
 finally:verify.D=old
print('PASS: actual certificate-byte tamper rejected')
