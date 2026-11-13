import sys
if sys.flags.optimize: raise RuntimeError('Run verification without -O')
"""Replay recipe, exact exchange and full geometry. No optimizer or network."""
from pathlib import Path
import numpy as np,json,hashlib,argparse,subprocess,sys
from verify_geometry import verify,require
ROOT=Path(__file__).resolve().parent;D=ROOT/'data'

def keys(a):return {r.tobytes() for r in np.asarray(a,dtype='<i2')}
def run(d):
    n,old=(231814,231776) if d==29 else (248180,248084)
    subprocess.run([sys.executable,str(ROOT/f'verify_recipe{d}.py')],check=True)
    f=np.load(D/f'coordinates{d}_{n}.npz',allow_pickle=False)['points'].astype(np.int16)
    source=np.load(D/f'source{d}_{old}.npz',allow_pickle=False)['points'].astype(np.int16)
    require(np.all(f%2==0) and np.all(source%2==0),'Integral norm32 scaling')
    f=f//2;source=source//2
    with np.load(D/f'delta{d}_{n}.npz',allow_pickle=False) as z:
        added=z['added'];removed=z['removed']
    sk,ak,rk,fk=keys(source),keys(added),keys(removed),keys(f)
    require(len(sk)==old and len(fk)==n,'Source/final distinct counts')
    require(len(ak)==len(added) and len(rk)==len(removed),'Distinct delta')
    require(rk<=sk and not ak&sk,'Delta membership')
    require((sk-rk)|ak==fk,'Exact source/delta reconstruction')
    E=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,d)]
    module=np.load(D/('oriented29_heads.npz' if d==29 else 'native30module.npz'))['points']
    require(keys(added[:,E])<=keys(module),'Every addition comes from proved module')
    require(np.all(added[:,np.setdiff1d(np.arange(d),E)]==0),'Module embedding')
    # Different arithmetic implementation: all changed pairs using direct int64 sums.
    retained=np.array([r for r in source if r.tobytes() not in rk],dtype=np.int64)
    for row in added.astype(np.int64):
        dot=np.zeros(len(retained),dtype=np.int64)
        for k in np.flatnonzero(row):dot+=row[k]*retained[:,k]
        require(np.max(dot)<=16,'Added versus retained')
    gram=added.astype(np.int64)@added.T.astype(np.int64);np.fill_diagonal(gram,-32)
    require(np.all(gram<=16),'Added pairs')
    result=verify(D/f'coordinates{d}_{n}.npz',precision=64)
    result.update(added=len(added),removed=len(removed),previous_count=old,exchange='PASS',recipe='PASS')
    print(json.dumps(result,indent=2));return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dimension',type=int,choices=[29,30]);args=ap.parse_args()
    manifest=json.loads((D/'manifest.json').read_text())
    for name,digest in manifest.items():require(hashlib.sha256((D/name).read_bytes()).hexdigest()==digest,'File hash '+name)
    for d in ([args.dimension] if args.dimension else [29,30]):run(d)
