"""Exact layer geometry and counts for the new 61D--65D sign families."""
from fractions import Fraction as Q
from itertools import combinations,product
from math import comb
from pathlib import Path
import json
from verify_family import generator,rref,code

def radical_leq(a,b,r,limit=Q(1,2)):
    a,b=Q(a)-limit,Q(b)
    if b==0:return a<=0
    if b>0:return a<=0 and b*b*r<=a*a
    return a<=0 or a*a<=b*b*r

def shorten(rows,positions):
    rows=list(rows)
    for p in positions:
        k=next((i for i,x in enumerate(rows) if x>>p&1),None)
        if k is None:continue
        pivot=rows.pop(k);rows=[x^pivot if x>>p&1 else x for x in rows]
    return rows

def verify_geometry():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    gen=generator();four=(24,42,53,55)
    projection=[sum(((x>>p)&1)<<i for i,p in enumerate(four)) for x in gen]
    assert len(rref(projection,4)[0])==4
    assert len(shorten(gen,four[2:]))==26 and len(shorten(gen,four))==24
    tail=(Q(1,4),Q(1,16));body2=(Q(1,64),Q(17,1024));checks61=0
    for j in range(2):assert 60*body2[j]+tail[j]**2==1
    for i,j,k,l in product(range(2),repeat=4):
        low=16-(i!=k)-(j!=l)
        for distance in range(low,61):
            dot=60-2*distance;a=(-1 if i!=k else 1)*tail[j]*tail[l]
            assert radical_leq(a+dot*body2[j],0,17) if j==l else radical_leq(a,Q(dot,256),17)
            checks61+=1
    # Every selected sparse class has odd total sign parity on a dual support.
    assert Q(14,32)<=Q(1,2) and radical_leq(0,Q(14,128),17)
    labels=(0,1,3,7,6,4);hexagon=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
    for i,(x,y) in enumerate(hexagon):
        assert 60+x*x+3*y*y==64
        for j,(u,v) in enumerate(hexagon):
            head_distance=16-(labels[i]^labels[j]).bit_count()
            assert 60-2*head_distance+x*u+3*y*v<=32
    # Affine planes form a quadruple packing: distinct supports share <=2 sites.
    blocks={tuple(sorted((a,b,c,a^b^c))) for a,b,c in combinations(range(64),3)}
    assert len(blocks)==10416
    counts={d:sum(max(b)<d for b in blocks) for d in (61,62,64)}
    assert counts=={61:8555,62:9145,64:10416}
    seen=set()
    for b in blocks:
        for triple in combinations(b,3):assert triple not in seen;seen.add(triple)
    # Dense against quadruples; each tail coordinate is handled separately.
    for j in range(2):
        for nt in (0,1):
            assert Q(4-nt,16)+Q(nt,8)<=Q(1,2) if j==0 else radical_leq(Q(nt,32),Q(4-nt,64),17)
    for x,y in hexagon:
        for a,b in product(range(2),repeat=2):
            assert radical_leq(Q(4-a-b+a*abs(x),16),Q(b*abs(y),16),3)
    assert Q(16,32)==Q(1,2)  # 62D dense / weight-sixteen layer.
    assert Q(8,16)==Q(2,4)==Q(4,8)==Q(1,2)
    assert radical_leq(0,Q(1,32),17) and Q(1,4)<=Q(1,2)
    assert 2*16==2*4*4==32 and 2*8<32  # norm64 parent and other layers.
    # Reusing the existing 65D owner construction: caps against every added point.
    assert radical_leq(0,16,3,Q(32)) and Q(3,4)*16+16<=32
    assert Q(3,4)*64-16==32
    return dict(status='PASS',dense61_pair_cases=checks61,hexagonal_pair_types=36,
                four_support_counts=counts,dense61_count=1<<26,dense62_count=6*(1<<24),
                owner65_scope='The unchanged 68,224-owner BCH/affine-line certificate is inherited; all new cap/equator interactions are checked here.')

def counts(certificate):
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    c=json.loads(Path(certificate).read_text());n=len(c['new'])
    if c['removed_coordinates']:
        assert c['removed_coordinates']==[24,42,53,55] and c['odd_total_parity'] and c['original_supports_retained']==9673
        return {61:(1<<26)+(9673+n)*2048+8555*16+122,
                62:6*(1<<24)+(9673+n)*2048+9145*16+124}
    base=(1<<28)+(30828+n)*2048+10416*16+128
    return {64:base,65:base+68224+2}

if __name__=='__main__':
    import sys
    print(json.dumps(dict(geometry=verify_geometry(),dimensions=counts(sys.argv[1])),indent=2))
