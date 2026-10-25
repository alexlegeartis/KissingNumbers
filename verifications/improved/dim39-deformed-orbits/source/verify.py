"""Exact complete orbit verifier. Requires Python and NumPy; run without -O.

Every geometric check uses integers. The Walsh transform evaluates every
relative group element, including pairs far from the original binary code.
The construction consists of exact normalized integer vectors, not rounded
unit vectors. No optimization output or cached distance table is trusted.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
import verify_base

HERE=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise ValueError(message)

def group_coordinates(group):
    require(all(type(w) is int and 0<=w<1<<39 for w in group),'39-bit group words')
    require(len(set(group))==len(group) and 0 in group,'Distinct group elements')
    pivots={};basis=[]
    for word in group:
        residual=word
        while residual:
            bit=residual.bit_length()-1
            if bit in pivots:residual^=pivots[bit]
            else:
                pivots[bit]=residual;basis.append(word);break
    enumeration=[0]
    for b in basis:enumeration += [w^b for w in enumeration]
    require(len(enumeration)==len(group) and set(enumeration)==set(group),'Exactly a binary linear group')
    characters=[sum(((b>>j)&1)<<k for k,b in enumerate(basis)) for j in range(39)]
    return basis,enumeration,characters

def orbit_products(V,first,second,characters,order):
    values=np.zeros((len(first),order),dtype=np.int64)
    for j,c in enumerate(characters):values[:,c]+=V[first,j]*V[second,j]
    stride=1
    while stride<order:
        blocks=values.reshape(len(first),-1,2*stride)
        left=blocks[:,:,:stride].copy();right=blocks[:,:,stride:].copy()
        blocks[:,:,:stride]=left+right;blocks[:,:,stride:]=left-right
        stride*=2
    return values

def main():
    if not __debug__:raise RuntimeError('Run without -O.')
    start=time.monotonic();verify_base.main()
    cert=json.loads((HERE/'orbit_certificate.json').read_text())
    require(cert['dimension']==39,'Dimension')
    require(cert['norm_lower']==99996 and cert['norm_upper']==100004,'Norm bounds')
    lo=cert['norm_lower'];hi=cert['norm_upper'];rows=cert['representatives']
    require(len(rows)>0 and all(len(v)==39 for v in rows),'Representative dimensions')
    require(all(type(x) is int and abs(x)<=100000 for v in rows for x in v),'Integer coefficients and safety bound')
    V=np.asarray(rows,dtype=np.int64);n2=np.sum(V*V,axis=1)
    require(np.all((lo*lo<=n2)&(n2<=hi*hi)),'Exact squared norms')
    # All Walsh partial sums are bounded by the sum of absolute coordinate
    # products, at most 39*100000^2. Sparse checks have sums <=8*100000.
    require(2*39*100000**2<2**63 and 2*(8*100000)**2<2**63,'Signed int64 arithmetic safety')
    mag=np.sort(abs(V),axis=1)[:,::-1]
    top8=mag[:,:8].sum(1);top2=mag[:,:2].sum(1);heavy=2*mag[:,0]+mag[:,1:5].sum(1)
    require(np.all(top8*top8<=2*lo*lo),'All signed eight-support vectors')
    require(np.all(2*top2*top2<=lo*lo),'All signed root vectors')
    require(np.all(heavy*heavy<=2*lo*lo),'All signed heavy-star vectors')
    basis,group,characters=group_coordinates(cert['group_words']);order=len(group)
    require(order<=65536,'Verifier resource bound')
    # Independent direct dot-product cross-check for every group action on
    # these deterministic pairs, including self-pairs and the inserted seed.
    signs=np.array([[1-2*((g>>j)&1) for j in range(39)] for g in group],dtype=np.int64)
    probes=sorted(set([(0,0),(0,len(V)-1),(len(V)-1,len(V)-1),(len(V)//2,len(V)-1)]))
    for a,b in probes:
        transformed=orbit_products(V,np.array([a]),np.array([b]),characters,order)[0]
        require(np.array_equal(transformed,signs@(V[a]*V[b])),'Direct/transform cross-check')
    first,second=np.triu_indices(len(V));maximum=-2**62;checked=0
    witness=None
    for begin in range(0,len(first),1024):
        a=first[begin:begin+1024];b=second[begin:begin+1024]
        dots=orbit_products(V,a,b,characters,order)
        # Identity is the only excluded pair: a point with itself.
        same=a==b;dots[same,0]=-2**62
        local=int(dots.max())
        if local>maximum:
            maximum=local;row,g=np.unravel_index(dots.argmax(),dots.shape)
            witness=[int(a[row]),int(b[row]),group[g]]
        require(local<=lo*lo//2,'A distinct orbit pair exceeds 1/2')
        checked+=len(a)*order-int(same.sum())
    base=json.loads((HERE/'base_verification.json').read_text())
    require(base['all_checks_passed'],'Sparse certificate verification')
    sparse=base['middle_vectors']+base['retained_roots']+base['added_heavy_vectors']
    require(sparse==441244,'Sparse layer size')
    total=sparse+len(V)*order
    require(total==cert['claimed_size'],'Total count')
    report={'all_checks_passed':True,'dimension':39,'lower_bound':total,
        'sparse_points':sparse,'orbit_representatives':len(V),'group_order':order,
        'dense_points':len(V)*order,'complete_relative_orbit_products_checked':checked,
        'maximum_integer_inner_product':maximum,'maximum_witness':witness,
        'allowed_integer_inner_product':lo*lo//2,'norm_lower':lo,'norm_upper':hi,
        'top8_square_slack':int((2*lo*lo-top8*top8).min()),
        'top2_square_slack':int((lo*lo-2*top2*top2).min()),
        'heavy_square_slack':int((2*lo*lo-heavy*heavy).min()),
        'arithmetic':'Bounded integer arithmetic; complete relative orbit-pair enumeration',
        'certificate_sha256':hashlib.sha256((HERE/'orbit_certificate.json').read_bytes()).hexdigest(),
        'elapsed_seconds':time.monotonic()-start}
    (HERE/'orbit_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));print('VERIFIED: tau_39 >= '+str(total)+'.')

if __name__=='__main__':main()
