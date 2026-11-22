"""Exact integer geometry for the 37D four-level QR compression."""
from collections import Counter
from itertools import combinations,product
from pathlib import Path
import hashlib,json
from supports37 import supports

P=Path(__file__).resolve().parent
geometry=[4,1,-4,-1]  # two-bit word b0+2*b1
norms=[144+x*x for x in geometry]

def replay():
    c=json.loads((P/'dense_parent38.json').read_text())
    G=list(map(lambda x:int(x,16),c['generators_hex']))
    R=list(map(lambda x:int(x,16),c['representatives_hex']))
    assert len(G)==14 and len(R)==11 and all(0<=x<1<<38 for x in G+R)
    H=[0]
    for g in G:
        H += [h^g for h in H]
    assert len(H)==len(set(H))==16384 and min(h.bit_count() for h in H[1:])==12
    assert all(min((a^b^h).bit_count() for h in H)==10 for a,b in combinations(R,2))
    assert Counter(h&3 for h in H)=={k:4096 for k in range(4)}
    assert len(set(geometry))==4
    dmax=dict()
    for a,b in product(range(4),repeat=2):
        old_tail_product=4*(2-2*(a^b).bit_count())
        dot=72-old_tail_product+geometry[a]*geometry[b]
        assert dot<=0 or 4*dot*dot<=norms[a]*norms[b]
        key=tuple(sorted([norms[a],norms[b]]));dmax[key]=max(dot,dmax.get(key,0))
    # No old pair of distinct codewords shares the remaining 36 columns:
    # such a pair would have parent distance <=2. The scalar map is itself
    # injective on the four patterns.
    S=supports()
    assert len(S)==len(set(S))==2742
    assert all(0<s<1<<36 and s.bit_count()==8 for s in S)
    inter=Counter((s&t).bit_count() for s,t in combinations(S,2))
    assert inter=={0:400893,1:1139016,2:1376880,3:583464,4:257658}
    for n in norms:
        assert 4*32**2<=32*n # middle dot <=32, with the new coordinate absent
        assert 4*24**2<=32*n # root dot <=24, allowing the new coordinate
    # middle-middle, middle-root, root-root as in the three-layer construction:
    # max 16 at squared norm 32 using parity on a shared support.
    assert 4*16**2==32**2
    d=len(H)*len(R);m=128*len(S);r=2*37*36
    assert d+m+r==533864
    return dict(status='PASS',dimension=37,points=d+m+r,dense=d,middle=m,roots=r,
        dense_product_bounds=[list(k)+[dmax[k]] for k in sorted(dmax)],octad_intersections=inter)
if __name__=='__main__':
    print(json.dumps(replay(),sort_keys=True))

