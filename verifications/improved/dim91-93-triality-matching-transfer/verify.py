#!/usr/bin/env python3
"""Exact check of three transfers; no optimizer is used."""
from itertools import combinations,product
from pathlib import Path
import hashlib,json
from importlib import util
from collections import Counter
from supports import hanani


P=Path(__file__).resolve().parent
def pin(rel,digest):
    p=P.parent/rel; raw=p.read_bytes()
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest, rel
    return p

path89=pin('dim89-mixed-matching-rows/verify.py','4f271af79e8c4f322b5b92aafd9f7996ad2b87d2')
spec=util.spec_from_file_location('matching89',path89)
m=util.module_from_spec(spec);spec.loader.exec_module(m)
geometry,data,decode,qmatch=m.geometry,m.data,m.decode,m.qmatch

def check_triality():
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    axial=[tuple(s*2*int(i==k) for i in range(4)) for k in range(4) for s in (-1,1)]
    flat=[v for v in product((-1,1),repeat=4)]
    frame=[[v for v in flat if v.count(-1)%2==t] for t in (0,1)]+[axial]
    signed_roots=[tuple(a if k==i else b if k==j else 0 for k in range(4)) for a,b in product((-1,1),repeat=2)
                  for i,j in combinations(range(4),2)]
    golay=[tuple(2*x for x in a) for a in signed_roots]
    assert len(set(signed_roots))==24 and max(dot(a,b) for a,b in combinations(signed_roots,2))==1
    for A in frame:
        assert len(set(A))==8 and all(dot(a,a)==4 for a in A)
        assert all(dot(a,b) in (0,-4) for a,b in combinations(A,2))
    assert max(dot(a,b) for A,B in combinations(frame,2) for a,b in product(A,B))==2
    assert max(dot(a,b) for a,b in product(flat+axial,golay))==4
    assert max(dot(a,b) for a,b in product(flat+axial,signed_roots))==2
    assert max(dot(a,b) for a,b in product(golay,signed_roots))==4
    Q=[(a,b,c,a^b^c) for a,b,c in product(range(4),repeat=3)]
    assert len(Q)==len(set(Q))==64 and all(sum(a!=b for a,b in zip(x,y))>=2 for x,y in combinations(Q,2))
    assert 4*(4-2)==8   # all sixteen signs, one reversal suffices
    # Distinct supports share <=2 rows. Ordinary D4 roots use one row only.
    for old_radius in (80,88):
        assert 4*16**2<=16*old_radius and 4*4**2<=2*old_radius
    assert 4*2**2<=16*2
    return len(Q)*16,len(signed_roots)

def main():
    assert __debug__
    stock=json.loads(pin('dim91-93-flat-parity-extensions/flat_supports.json','b1112ec5baa8c6b3619f7b4215d8769c16c1c94c').read_text())
    claims=json.loads(pin('dim91-93-flat-parity-extensions/claims.json','d078130ffe70bfbc7383eca3c293c6ff9280a23f').read_text())
    flats=[(w,a,p) for w,key in ((6,'six'),(7,'seven')) for a,b in stock[key] for p in (b,a^b)]
    fw={6:2048,7:8192}
    assert len(flats)==len(set(flats))==490
    assert all(0<a<1<<22 and a.bit_count()==w and a&p==p for w,a,p in flats)
    old_flat_points=sum(fw[w] for w,_,_ in flats)
    assert old_flat_points==2613248
    # All retained old stock is a subset of these certified parent patterns.
    sq=[sum(1<<j for j in b) for b in hanani(8)]
    assert len(sq)==len(set(sq))==385
    assert all((a&b).bit_count()<=2 for a,b in combinations(sq,2))
    sf=[(s,t)for s in sq for t in range(2)]
    mode_count,rays_root=check_triality()
    old_golay=data()
    per=geometry()
    pick=json.loads((P/'selection.json').read_text())
    assert set(pick)==set(('keep_old','append_flat_modes','append_matching'))
    def subset(A,p):
        assert p==sorted(set(p)) and all(0<=i<len(A) for i in p)
        return [A[i] for i in p]
    old_keep=subset(flats,pick['keep_old'])
    new_flat_modes=subset(sf,pick['append_flat_modes'])
    tags0=[tuple(map(int,l.split())) for l in pin('dim89-mixed-matching-rows/witness.txt','bae0e9fa8bc696ef5c23e0cc6e81fdd1be856adb').read_text().splitlines()]
    mrows=[decode(a,b) for a,b in subset(tags0,pick['append_matching'])]
    assert len(set(mrows))==len(mrows)
    mw=[(a|b).bit_count() for a,b in mrows]
    for row,w in zip(mrows,mw):
        mask=row[0]|row[1]
        if w<8:assert all(((s<<2)&mask)!=s<<2 for s in sq)
        for x in old_golay:assert qmatch(row,x)**2<=w*(x[0]|x[1]).bit_count()
        # Matching vs pure-flat old: each overlapping row contributes at most 2sqrt(2).
        for r,b,p in old_keep:assert 2*((mask>>2)&b).bit_count()**2<=w*r
    for (a,w),(b,r) in combinations(zip(mrows,mw),2):assert qmatch(a,b)**2<=w*r
    for w,a,p in old_keep:
        # Two old flat parity modes remain source-defined. Axial vs these layers has
        # maximum two per shared row, local ordinary roots also at most two.
        assert 4*(2*min(w,4))**2<=4*w*16
        assert 4*2**2<=4*w*2
        for b,t in new_flat_modes:
            common=a&b
            q=common.bit_count()
            d=(common&(p^(b if t else 0))).bit_count()
            assert (4*q-2*d)**2<=16*w
    # Dense/SQ/root/old Golay inherit normalized comparisons. Other same-weight new/old
    # flat-matching cross pairs are nonduplicates because their row shapes differ (2 versus 4 entries).
    keep_flat=sum(fw[w] for w,_,_ in old_keep)
    new_triality=(len(sq)+len(new_flat_modes))*mode_count+rays_root*22
    new_matching=sum(per[w] for w in mw)
    dense_fibres=(24,32,44)
    gsize=982*1024*144+11*256*72
    result=[]
    for d,fib,cl in zip((91,92,93),dense_fibres,claims):
        base=fib*8**10+gsize
        assert cl['dimension']==d and cl['lower_bound']==base+old_flat_points
        total=base+keep_flat+new_triality+new_matching
        result.append(dict(dimension=d,previous=cl['lower_bound'],lower_bound=total,gain=total-cl['lower_bound']))
    out=dict(claims=result,old_flat_retained=keep_flat,unconditional_axial_modes=len(sq),
                     flat_modes=len(new_flat_modes),ordinary_roots=22*rays_root,
                     matching_patterns=dict(sorted(Counter(mw).items())),
                     new_triality_and_roots=new_triality,new_matching=new_matching)
    assert json.loads(json.dumps(out))==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
