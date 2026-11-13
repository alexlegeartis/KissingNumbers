from pathlib import Path
import numpy as np
from verify_geometry import verify
if not __debug__: raise RuntimeError('Run without -O: exact certificate assertions are required.')
from check_module_recipe import main as check_module
check_module()
P=Path(__file__).resolve().parent/'data'
x=np.load(P/'source251254.npz')['points'];f=np.load(P/'coordinates31_251420.npz')['points']//2
d=np.load(P/'delta251420.npz');add=d['added'];rem=d['removed']
assert len(x)==251254 and len(f)==251420 and add.shape==(2218,31) and rem.shape==(2052,31)
xk={r.tobytes() for r in x};rk={r.tobytes() for r in rem};assert rk<=xk
r=np.r_[np.array([a for a in x if a.tobytes() not in rk]),add]
assert {a.tobytes() for a in r}=={a.tobytes() for a in f} and len(r)==len(f)
# Independent sparse integer exact changed pairs.
for row in add:
 products=np.zeros(len(f),np.int32)
 for k in np.flatnonzero(row):products+=int(row[k])*f[:,k].astype(np.int32)
 assert np.sum(products==32)==1 and np.max(products[products<32])<=16
print('Exact int32 exchange replay PASS')
print(verify(P/'coordinates31_251420.npz',precision=64))
