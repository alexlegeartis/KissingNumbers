"""Portable exact verifier for the mixed 31D certificates; NumPy only.

For midpoint rows, physical coordinates are (x_0,...,x_29,sqrt(2)*x_30).
For centre rows, the first 24 coordinates are unchanged; their tail uses
the same metric followed by multiplication by sqrt(q), where q is the
certificate's positive rational centre_scale_squared (default 2).
The dual-centre certificates use q=8/3. Pure tail rows use the midpoint
tail metric. All points have norm squared 32 and maximum product 16.

Integer dot products use float64 BLAS exactly: small integer operands and
bounded absolute partial sums, with no rounding tolerance. Radical products
are compared by integer squaring, including their signs.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np


def leech_shell():
    generator=np.array([1,0,1,0,1,1,1,0,0,0,1,1],dtype=np.uint8)
    matrix=np.zeros((12,23),dtype=np.uint8)
    for i in range(12):matrix[i,i:i+12]=generator
    bits=((np.arange(4096)[:,None]>>np.arange(12))&1).astype(np.uint8)
    golay=np.zeros((4096,24),dtype=np.uint8)
    golay[:,:23]=(bits@matrix)%2;golay[:,23]=golay[:,:23].sum(1)%2
    assert dict(zip(*np.unique(golay.sum(1),return_counts=True)))=={0:1,8:759,12:2576,16:759,24:1}
    rows=[]
    for i,j in itertools.combinations(range(24),2):
        for a,b in itertools.product([4,-4],repeat=2):
            row=np.zeros(24,dtype=np.int8);row[i]=a;row[j]=b;rows.append(row)
    signs=np.array([[(-2 if k>>i&1 else 2) for i in range(8)] for k in range(256) if bin(k).count('1')%2==0],dtype=np.int8)
    for word in golay[golay.sum(1)==8]:
        block=np.zeros((128,24),dtype=np.int8);block[:,np.flatnonzero(word)]=signs;rows.extend(block)
    sign_words=np.where(golay,-1,1).astype(np.int8)
    for j in range(24):
        block=sign_words.copy();block[:,j]*=-3;rows.extend(block)
    m=np.array(rows,dtype=np.int8)
    assert m.shape==(196560,24) and np.all((m.astype(int)**2).sum(1)==32)
    assert len({row.tobytes() for row in m})==len(m)
    return m


def integer_dot(a,b):
    result=a.astype(np.float64)@b.astype(np.float64).T
    assert np.all(result==np.rint(result))
    return result.astype(np.int32)


def radical_leq(a,b,threshold=16,num=2,den=1):
    """a + sqrt(num/den)*b <= threshold, exactly for integer a,b."""
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64);limit=threshold-a
    return np.where(b>=0,(limit>=0)&(num*b*b<=den*limit*limit),(limit>=0)|(num*b*b>=den*limit*limit))


def verify(path):
    start=time.time();z=np.load(path,allow_pickle=False)
    for name in z.files:assert np.issubdtype(z[name].dtype,np.integer)
    x=z['midpoints'].astype(np.int32);c=z['centres'].astype(np.int32);pure=z['pure_tails'].astype(np.int32)
    removed=z['removed_leech'].astype(np.int32);ends=z['midpoint_owner_indices'].astype(np.int64)
    ce=z['centre_owner_pairs'].astype(np.int32)
    num,den=map(int,z['centre_scale_squared']) if 'centre_scale_squared' in z else (2,1)
    assert (num,den) in [(2,1),(8,3)]
    assert x.ndim==2 and x.shape[1]==31 and c.ndim==2 and c.shape[1]==31 and pure.shape[1]==7
    w=np.array([1]*6+[2],dtype=np.int32);fullw=np.r_[np.ones(24,dtype=np.int32),w]
    m=leech_shell();keys={row.tobytes() for row in m}
    assert np.all(abs(removed)<=4) and len(np.unique(removed,axis=0))==len(removed)
    assert all(row.astype(np.int8).tobytes() in keys for row in removed)
    assert len(np.unique(x,axis=0))==len(x) and len(np.unique(c,axis=0))==len(c)
    assert len(np.unique(pure,axis=0))==len(pure)
    assert np.all((x*x*fullw).sum(1)==32) and np.all((x[:,24:]**2*w).sum(1)==8)
    assert ends.shape==(len(x),2) and np.all(ends>=0) and np.all(ends<len(removed))
    assert np.all(ends[:,0]!=ends[:,1]) and np.all(2*x[:,:24]==removed[ends].sum(1))
    assert np.all((c[:,:24]**2).sum(1)==16) and np.all(num*(c[:,24:]**2*w).sum(1)==16*den)
    assert ce.shape==(len(c),2,24) and np.all(abs(ce)<=4)
    assert all(row.astype(np.int8).tobytes() in keys for row in ce.reshape(-1,24))
    assert np.all(ce.sum(1)==2*c[:,:24]) and np.all((ce[:,0]*ce[:,1]).sum(1)==0)
    # Every centre head clears all Leech points, including its witnesses:
    # <(u+v)/2,u>=(32+0)/2=16. Other points clear by convexity.
    assert np.all((pure*pure*w).sum(1)==32)
    pg=integer_dot(pure*w,pure);np.fill_diagonal(pg,-100);assert pg.max(initial=-100)<=16
    cg=den*integer_dot(c[:,:24],c[:,:24])+num*integer_dot(c[:,24:]*w,c[:,24:])
    np.fill_diagonal(cg,-100);assert cg.max(initial=-100)<=16*den
    cp=integer_dot(c[:,24:]*w,pure);assert radical_leq(0,cp,num=num,den=den).all()
    max_mid=-100
    for at in range(0,len(x),384):
        block=x[at:at+384];last=at+len(block)
        dot=integer_dot(block*fullw,x[:last]);dot[np.arange(len(block)),at+np.arange(len(block))]=-100
        assert dot.max()<=16,(at,int(dot.max()));max_mid=max(max_mid,int(dot.max()))
        assert integer_dot(block[:,24:]*w,pure).max(initial=-100)<=16
        aa=integer_dot(block[:,:24],c[:,:24]);bb=integer_dot(block[:,24:]*w,c[:,24:])
        assert radical_leq(aa,bb,num=num,den=den).all(),('midpoint-centre',at)
        if at%15360==0:print('checked midpoint rows',last,'of',len(x),flush=True)
    total=196560-len(removed)+len(x)+len(c)+len(pure)
    result=dict(status='PASS EXACT',dimension=31,points=total,removed=len(removed),inserted_midpoints=len(x),inserted_centres=len(c),pure=len(pure),
                maximum_midpoint_product=max_mid,centre_scale_squared=[num,den],arithmetic='Exact small-integer BLAS and signed integer squaring for quadratic radicals; no tolerance.',
                certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),seconds=time.time()-start)
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);ap.add_argument('--output',type=Path);args=ap.parse_args()
    result=verify(args.certificate)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
