"""Exact root-star extensions; standard library only, no search is trusted."""
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent


def verify(path):
    start=time.monotonic();trade=json.loads(path.read_text())
    assert trade['format']=='root_star_trade_v1'
    n=trade['dimension'];assert type(n) is int and n in (33,34,35)
    expected='source/multi35/certificate35.json' if n==35 else 'source/classic/certificate%d.json'%n
    assert trade['baseline']==expected
    base=PACKAGE/expected
    assert hashlib.sha256(base.read_bytes()).hexdigest()==trade['baseline_sha256']
    spec=importlib.util.spec_from_file_location('root_trade_baseline',base.parent/'verify.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    original=module.verify(base)
    assert original.get('all_checks_pass',original.get('verified',False))
    c=json.loads(base.read_text());assert c['dimension']==n
    supports=c.get('full_parity_supports',[(s,0) for s in c.get('full_even_supports',[])])
    exceptional=c.get('exceptional_antipodal_lines',[])
    old_removed=set(map(tuple,c.get('removed_root_lines',[])))
    heads=c.get('head_antipodal_lines',[])
    coordinates=c.get('dense_coordinates',list(range(32)))
    if n==34:
        a,b=c['heavy_pair_extension']['coordinates'];old_removed={(a,b,-1)}
        heads=[]
        for sign in (-1,1):
            q=[1]*n;q[a],q[b]=3*sign,-3*sign;heads.append(q)
    required=set();added=[]
    for star in trade['stars']:
        a=star['center'];leaves=star['leaves']
        assert type(a) is int and 0<=a<n
        assert len(leaves)==len(set(leaves))==4 and a not in leaves
        assert all(type(b) is int and 0<=b<n for b in leaves)
        for b in leaves:
            i,j=sorted((a,b))
            required.update((i,j,sign) for sign in (-1,1))
        for center_sign in (-1,1):
            for mask in range(16):
                if bin(mask).count('1')%2:continue
                v=[0]*n;v[a]=4*center_sign
                for k,b in enumerate(leaves):v[b]=2-4*((mask>>k)&1)
                added.append(tuple(v))
    assert len(added)==len(set(added))==16*len(trade['stars'])
    assert all(sum(x*x for x in v)==32 for v in added)
    removed=set(map(tuple,trade['removed_root_lines']))
    assert len(removed)==len(trade['removed_root_lines']) and removed==required-old_removed
    roots=[(a,b,sign) for a,b in itertools.combinations(range(n),2)
           for sign in (-1,1) if (a,b,sign) not in removed|old_removed]
    maxima=dict(middle=0,exceptional=0,roots=0,stars=0,heads=0,dense_bound=0)
    for v in added:
        dense_bound=sum(abs(v[j]) for j in coordinates)
        assert dense_bound<=16
        maxima['dense_bound']=max(maxima['dense_bound'],dense_bound)
        for mask,parity in supports:
            values=[v[j] for j in range(n) if mask>>j&1]
            score=sum(abs(x) for x in values)
            if sum(x<0 for x in values)%2!=parity:score-=2*min(map(abs,values))
            score*=2
            assert score<=16
            maxima['middle']=max(maxima['middle'],score)
        for mask,minus in exceptional:
            dot=abs(2*sum(v[j]*(1-2*((minus>>j)&1)) for j in range(n) if mask>>j&1))
            assert dot<=16
            maxima['exceptional']=max(maxima['exceptional'],dot)
        for a,b,sign in roots:
            dot=abs(4*(v[a]+sign*v[b]));assert dot<=16
            maxima['roots']=max(maxima['roots'],dot)
        for q in heads:
            dot=abs(sum(x*y for x,y in zip(v,q)))
            assert 4*dot*dot<=32*sum(x*x for x in q)
            maxima['heads']=max(maxima['heads'],dot)
    for v,w in itertools.combinations(added,2):
        dot=sum(x*y for x,y in zip(v,w));assert dot<=16
        maxima['stars']=max(maxima['stars'],dot)
    total=c['points']-2*len(removed)+len(added)
    assert total==trade['points'] and total>c['points']
    return dict(dimension=n,points=total,baseline_points=c['points'],
                stars=len(trade['stars']),new_star_points=len(added),removed_root_points=2*len(removed),
                maximum_integer_products=maxima,all_checks_pass=True,
                certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                arithmetic='Python integers; exact normalized head comparisons',
                elapsed_seconds=time.monotonic()-start)


if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O.')
    records=[verify(path) for path in sorted(HERE.glob('certificate*.json'))]
    (HERE/'verification.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records,indent=2))
