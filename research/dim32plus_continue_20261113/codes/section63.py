"""Recount an exact 63D hyperplane section of an explicit 64D family."""
from itertools import combinations
from math import comb
from pathlib import Path
import json,sys
import numpy as np
from verify_family import SOURCE,code,generator,rref,span

def kernel(checks):
    rows,piv=rref(checks,16);free=[j for j in range(16) if j not in piv]
    return [(1<<j)^sum(((r>>j)&1)<<p for r,p in zip(rows,piv)) for j in free]

def local(s,h):
    return sum(((h>>c)&1)<<j for j,c in enumerate(c for c in range(64) if s>>c&1))

def verify_section(path,enumerate_words=True):
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    cert=json.loads(path.read_text());assert not cert['removed_coordinates']
    data={r['support']:code(r,False) for r in cert['old']+cert['new']}
    with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False)as z:old=list(map(int,z['words']))
    added=[r['support'] for r in cert['new']]
    octads=json.loads((SOURCE/'shared_parity/dual_octads.json').read_text());gen=generator()
    assert len(octads)==len(set(octads))==72
    def count(s,h):
        m=(s&h).bit_count()
        if m%2:return 0
        if m<8:return comb(m,m//2)*(1<<(11-m))
        assert m==8
        if s not in data:return 1120
        row=data[s]
        if h in row['dual']:return 1120 if not (row['phase']&h).bit_count()%2 else 0
        return 560
    candidates=[]
    for h in octads:
        assert h.bit_count()==8 and all((h&g).bit_count()%2==0 for g in gen)
        assert len(rref([local(h,g) for g in gen],8)[0])==7
        terms=dict(dense=70*(1<<21),old_sparse=sum(count(s,h) for s in old),new_sparse=sum(count(s,h) for s in added),four_supports=102900,poles=112)
        candidates.append(dict(octad=h,total=sum(terms.values()),terms=terms))
    best=max(candidates,key=lambda r:(r['total'],-r['octad']));h=best['octad']
    checks=0
    if enumerate_words:
        standard=[65535]+[sum(((i>>b)&1)<<i for i in range(16)) for b in range(4)]
        standard_words=span(kernel(standard));assert len(standard_words)==2048
        default_counts={m:sum((w&((1<<m)-1)).bit_count()==m//2 for w in standard_words) for m in (0,2,4,6,8)}
        by_support={r['support']:r for r in cert['old']+cert['new']}
        for s in old+added:
            m=(s&h).bit_count()
            if m%2:actual=0
            elif s not in by_support:actual=default_counts[m]
            else:
                r=by_support[s];words=span(kernel(r['checks']));assert len(words)==2048
                mask=local(s,h);actual=sum(((w^r['shift'])&mask).bit_count()==m//2 for w in words)
            assert actual==count(s,h);checks+=1
    # Choose the four-support labels independently, sending this octad to 0..7.
    blocks={tuple(sorted((a,b,c,a^b^c))) for a,b,c in combinations(range(64),3)}
    assert len(blocks)==10416
    four=sum(comb(m,m//2)*(1<<(4-m)) for block in blocks if (m:=sum(c<8 for c in block))%2==0)
    assert four==102900
    return dict(status='PASS',**best,octads_checked=72,sign_classes_independently_enumerated=checks,previous=184286468,gain=best['total']-184286468)

if __name__=='__main__':
    print(json.dumps(verify_section(Path(sys.argv[1]),'--no-enumerate' not in sys.argv[2:]),indent=2))
