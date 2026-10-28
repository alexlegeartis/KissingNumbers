"""Independent complete verification of an explicit integer kissing code.

NumPy is the only dependency. No baseline validity, search log, geometric
lemma, optimizer, or network connection is trusted. Every distinct pair
of antipodal lines is tested, with no floating-point tolerance.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def verify(path,block_size=128):
    start=time.time();raw=path.read_bytes()
    with np.load(path,allow_pickle=False) as z:
        assert int(z['squared_norm'])==32
        assert int(z['inner_product_limit'])==16
        assert int(z['coordinate_denominator'])==1
        a=z['points'];assert a.dtype.kind in 'iu' and a.ndim==2 and a.shape[1]==31
        assert np.all(a>=-8) and np.all(a<=8)
        points=a.astype(np.int16)
    assert np.all(np.sum(points.astype(np.int32)**2,axis=1)==32)
    unique=np.unique(points,axis=0);assert len(unique)==len(points),'Duplicate points'
    first=np.argmax(unique!=0,axis=1)
    positive=unique[np.arange(len(unique)),first]>0
    half=unique[positive];opposite=np.unique(-unique[~positive],axis=0)
    assert len(half)*2==len(points)
    assert np.array_equal(half,opposite),'Configuration is not antipodal'
    x=half.astype(np.float64);maximum=0;pairs=0;contacts=0
    for at in range(0,len(x),block_size):
        block=x[at:at+block_size];gram=block@x[at:].T
        rows,cols=np.tril_indices(len(block));gram[rows,cols]=0
        absolute=np.abs(gram);current=int(absolute.max())
        assert current<=16,('Forbidden inner product',at,current)
        maximum=max(maximum,current);contacts+=int(np.count_nonzero(absolute==16))
        pairs+=len(block)*(len(x)-at)-len(block)*(len(block)+1)//2
        if at%(64*block_size)==0:
            print('checked line rows',at,'/',len(x),'seconds',round(time.time()-start,1),flush=True)
    assert pairs==len(x)*(len(x)-1)//2
    # Each coordinate is an integer of magnitude <=8, so each product is an
    # integer <=64 in magnitude and each partial sum is <=31*64. Float64
    # represents every such product and sum exactly; reordering/FMA does not
    # alter the result. Thus the comparisons above are exact integer tests.
    return dict(status='PASS',dimension=31,points=len(points),squared_norm=32,
        inner_product_limit=16,antipodal_lines=len(half),distinct_line_pairs_checked=pairs,
        max_absolute_inner_product=maximum,contact_line_pairs=contacts,
        coordinate_sha256=hashlib.sha256(points.astype(np.int8).tobytes()).hexdigest(),
        file_sha256=hashlib.sha256(raw).hexdigest(),seconds=time.time()-start,
        arithmetic='Exact bounded integer products and sums in float64 BLAS; no tolerance.',
        scope='All norms, distinctness, antipodal partners, and all distinct line-pair inequalities; no external baseline assumptions.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coordinates',type=Path,nargs='?',default=ROOT/'data/coordinates31_251110.npz')
    p.add_argument('--report',type=Path);p.add_argument('--block-size',type=int,default=128)
    args=p.parse_args();assert args.block_size>0
    result=verify(args.coordinates,args.block_size)
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
