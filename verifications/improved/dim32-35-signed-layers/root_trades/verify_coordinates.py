"""Independent expanded-coordinate verification of every new-star pair type."""
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent


def verify(path):
    start=time.monotonic();t=json.loads(path.read_text());n=t['dimension']
    assert t['format']=='root_star_trade_v1' and n in (33,34,35)
    baseline='source/multi35/certificate35.json' if n==35 else 'source/classic/certificate%d.json'%n
    assert t['baseline']==baseline
    base=PACKAGE/baseline
    assert hashlib.sha256(base.read_bytes()).hexdigest()==t['baseline_sha256']
    spec=importlib.util.spec_from_file_location('coordinate_trade_baseline',base.parent/'verify_coordinates.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    receipt=module.verify(base)
    assert receipt.get('all_checks_pass',receipt.get('verified',False))
    c=json.loads(base.read_text());assert c['dimension']==n
    signed_supports=c.get('full_parity_supports',[(s,0) for s in c.get('full_even_supports',[])])
    patterns=np.array(list(itertools.product((-2,2),repeat=8)),dtype=np.int32)
    parity=np.sum(patterns<0,axis=1)%2
    layers=[]
    for mask,p in signed_supports:
        row=np.zeros((128,n),dtype=np.int32);indices=[j for j in range(n) if mask>>j&1]
        assert len(indices)==8 and p in (0,1)
        row[:,indices]=patterns[parity==p];layers.append(row)
    for mask,minus in c.get('exceptional_antipodal_lines',[]):
        v=np.array([2*(1-2*((minus>>j)&1)) if mask>>j&1 else 0 for j in range(n)],dtype=np.int32)
        layers.append(np.stack((v,-v)))
    middle=np.concatenate(layers)
    old_removed=set(map(tuple,c.get('removed_root_lines',[])))
    head_rows=c.get('head_antipodal_lines',[])
    if n==34:
        a,b=c['heavy_pair_extension']['coordinates'];old_removed={(a,b,-1)};head_rows=[]
        for sign in (-1,1):
            q=[1]*n;q[a],q[b]=3*sign,-3*sign;head_rows.append(q)
    required=set();stars=[]
    for entry in t['stars']:
        a=entry['center'];leaves=entry['leaves']
        assert type(a) is int and 0<=a<n and len(leaves)==len(set(leaves))==4 and a not in leaves
        assert all(type(j) is int and 0<=j<n for j in leaves)
        for b in leaves:
            i,j=sorted((a,b));required.update((i,j,sign) for sign in (-1,1))
        for center in (-4,4):
            for signs in itertools.product((-2,2),repeat=4):
                if sum(x<0 for x in signs)%2:continue
                v=[0]*n;v[a]=center
                for b,x in zip(leaves,signs):v[b]=x
                stars.append(v)
    H=np.array(stars,dtype=np.int32)
    assert len(H)==16*len(t['stars'])==len(np.unique(H,axis=0))
    assert np.all(np.sum(H*H,axis=1)==32)
    removed=set(map(tuple,t['removed_root_lines']))
    assert len(removed)==len(t['removed_root_lines']) and removed==required-old_removed
    roots=[]
    for a,b in itertools.combinations(range(n),2):
        for sa,sb in itertools.product((-4,4),repeat=2):
            if (a,b,sa*sb//16) in removed|old_removed:continue
            v=[0]*n;v[a],v[b]=sa,sb;roots.append(v)
    roots=np.array(roots,dtype=np.int32)
    words=np.array([0],dtype=np.uint32)
    for generator in c['dense_generator']:words=np.concatenate((words,words^np.uint32(generator)))
    dense=np.zeros((len(words),n),dtype=np.int32)
    coordinates=c.get('dense_coordinates',list(range(32)))
    dense[:,coordinates]=1-2*(((words[:,None]^np.uint32(c.get('dense_coset_shift',0)))
                               >>np.arange(32,dtype=np.uint32))&1).astype(np.int32)
    assert len(dense)==131072 and np.all(np.sum(dense*dense,axis=1)==32)
    maxima={};products=0
    for name,layer in [('dense',dense),('middle',middle),('roots',roots)]:
        dots=H@layer.T
        assert dots.max()<=16
        maxima[name]=int(dots.max());products+=dots.size
    gram=H@H.T;np.fill_diagonal(gram,0);assert gram.max()<=16
    maxima['stars']=int(gram.max())
    heads=np.array(head_rows,dtype=np.int32).reshape(-1,n)
    if len(heads):
        norms=np.sum(heads.astype(np.int64)**2,axis=1);dots=np.abs(H@heads.T).astype(np.int64)
        assert np.all(4*dots*dots<=32*norms)
        maxima['heads']=int(dots.max());products+=dots.size*2
    total=len(dense)+len(middle)+len(roots)+len(heads)*2+len(H)
    assert total==t['points'] and total==c['points']-2*len(removed)+len(H)
    return dict(dimension=n,points=total,baseline_points=c['points'],
                expanded_star_points=len(H),expanded_base_points=total-len(H),
                new_star_base_products=products,maximum_integer_products=maxima,
                all_checks_pass=True,certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                arithmetic='int32 coordinate products; int64 mixed-norm squares',
                elapsed_seconds=time.monotonic()-start)


if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O.')
    records=[verify(path) for path in sorted(HERE.glob('certificate*.json'))]
    (HERE/'coordinate_verification.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records,indent=2))
