"""Individual-ray bridge trades, with exact conflicts in Q(sqrt(2),sqrt(3)).

No bridge, axis, sign orbit, or six-dimensional block cardinality is fixed.
This is a finite candidate pool, not a parametrization of all spherical codes.
"""
from common import *
import itertools
import math
import time
import argparse
sys.path.insert(0,str(AUDIT))
from audit import clebsch_coefficients
from alternative_bridges import coefficients as alternative_coefficients

def sign_field(v):
    """Exact sign of a+b sqrt2+c sqrt3+d sqrt6, using rational enclosures."""
    v=[int(a) for a in v]
    if not any(v):return 0
    digits=16
    while True:
        scale=10**digits;lo=hi=v[0]*scale
        for coefficient,radicand in zip(v[1:],[2,3,6]):
            lower=math.isqrt(radicand*scale*scale);upper=lower+1
            lo+=coefficient*(lower if coefficient>=0 else upper)
            hi+=coefficient*(upper if coefficient>=0 else lower)
        if lo>0:return 1
        if hi<0:return -1
        digits*=2

def conflict(p,q,r,s,norm_product):
    if sign_field([p,q,r,s])<=0:return False
    square=[p*p+2*q*q+3*r*r+6*s*s,2*p*q+6*r*s,2*p*r+4*q*s,2*p*s+2*q*r]
    square=[4*v for v in square];square[0]-=norm_product
    return sign_field(square)>0

def make_pool(extended=False):
    rows=[];labels=[];lookup={}
    def add(c,label):
        c=np.asarray(c,dtype=np.int16);div=math.gcd(*[int(v) for v in c.ravel()]);assert div
        c=c//div;key=c.tobytes()
        if key not in lookup:lookup[key]=len(rows);rows.append(c);labels.append(label)
        return lookup[key]
    # Include both exact 840 configurations as positive controls.
    seeds=[]
    for name,func in [('standard',clebsch_coefficients),('alternative',alternative_coefficients)]:
        aa,bb=func();seed=[]
        for a,b in zip(aa,bb):
            c=np.zeros((4,12),dtype=np.int16);c[0]=a;c[1]=b;seed.append(add(c,name))
        seeds.append(seed)
    for block in [0,6]:
        for support_size in ([1,2,3,4,5,6] if extended else [1,2,4,6]):
            for support in itertools.combinations(range(6),support_size):
                for signs in itertools.product([-1,1],repeat=support_size):
                    c=np.zeros((4,12),dtype=np.int16);c[0,np.array(support)+block]=signs
                    add(c,f'block{block}_integer{support_size}')
        for special in range(6):
            for signs in itertools.product([-1,1],repeat=6):
                c=np.zeros((4,12),dtype=np.int16);c[0,block:block+6]=signs;c[0,block+special]*=2
                add(c,f'block{block}_norm9')
                # E6 spinor candidates: all signs offered individually. The
                # compatibility graph decides which parity classes can coexist.
                c=np.zeros((4,12),dtype=np.int16);c[0,block:block+6]=signs
                c[0,block+special]=0;c[2,block+special]=signs[special]
                add(c,f'block{block}_spinor')
        if extended:
            for shape in [(2,1,0,0,0,0),(2,1,1,1,1,0),(2,2,1,1,0,0)]:
                for absolute in sorted(set(itertools.permutations(shape))):
                    support=np.flatnonzero(absolute)
                    for signs in itertools.product([-1,1],repeat=len(support)):
                        c=np.zeros((4,12),dtype=np.int16)
                        c[0,block+support]=np.array(absolute)[support]*np.array(signs)
                        add(c,f'block{block}_integer_shape_{shape}')
    return np.array(rows),labels,np.array(seeds,dtype=np.uint16)

def graph(coeff):
    # Products are computed in binary64, but every operand/result is a small
    # integer, so these dot products are exact; integer conversion is asserted.
    a,b,c,d=[np.asarray(coeff[:,k],dtype=np.float64) for k in range(4)]
    products=[a@a.T+2*b@b.T+3*c@c.T+6*d@d.T,
              a@b.T+b@a.T+3*(c@d.T+d@c.T),
              a@c.T+c@a.T+2*(b@d.T+d@b.T),
              a@d.T+d@a.T+b@c.T+c@b.T]
    assert all(np.array_equal(v,np.rint(v)) for v in products)
    p,q,r,s=[v.astype(np.int16) for v in products];del products
    norm=np.diag(p).astype(np.int64)
    assert np.all(norm>0) and all(np.all(np.diag(v)==0) for v in [q,r,s])
    assert all(np.max(abs(v))<64 for v in [p,q,r,s])
    product=norm[:,None]*norm[None,:];assert product.max()<1024
    code=product.astype(np.uint64)
    for k,v in enumerate([p,q,r,s]):code|=(v.astype(np.int64)+64).astype(np.uint64)<<(10+7*k)
    unique,inverse=np.unique(code,return_inverse=True);table=[]
    for packed in unique:
        packed=int(packed);npair=packed&1023
        coeffs=[((packed>>(10+7*k))&127)-64 for k in range(4)]
        table.append(conflict(*coeffs,npair))
    bad=np.array(table,dtype=np.uint8)[inverse].reshape(len(norm),len(norm));np.fill_diagonal(bad,0)
    assert np.array_equal(bad,bad.T)
    return bad,norm,len(unique)

def main():
    p=argparse.ArgumentParser();p.add_argument('--extended',action='store_true');a=p.parse_args()
    reserve();name='trade_pool_extended' if a.extended else 'trade_pool'
    out=ROOT/'runs'/name;out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    coeff,labels,seeds=make_pool(a.extended);bad,norm,types=graph(coeff)
    assert all(len(set(seed))==840 and not bad[np.ix_(seed,seed)].any() for seed in seeds)
    raw=coeff[:,0]+np.sqrt(2)*coeff[:,1]+np.sqrt(3)*coeff[:,2]+np.sqrt(6)*coeff[:,3]
    x=raw/np.sqrt(norm[:,None]);g=x@x.T;check=abs(g-.5)>1e-12;np.fill_diagonal(check,False)
    assert np.array_equal(bad[check].astype(bool),(g>.5)[check])
    np.savez(out/'rays.npz',coeff=coeff,norm=norm,coordinates=x,seeds=seeds)
    np.array([len(coeff),len(seeds),840,842],dtype=np.uint32).tofile(out/'header.bin')
    bad.tofile(out/'weights.bin');seeds.tofile(out/'seeds.bin');save_json(out/'labels.json',labels)
    result=dict(vertices=len(coeff),edges=int(bad.sum()//2),distinct_exact_dot_tests=types,
                seed_sizes=[840,840],exact_seed_checks=True,
                elapsed_seconds=time.monotonic()-start,source_sha256=sha(__file__))
    save_json(ROOT/(name+'.json'),result);print(json.dumps(result),flush=True)

if __name__=='__main__':main()
