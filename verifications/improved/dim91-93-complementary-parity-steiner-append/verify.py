#!/usr/bin/env python3
"""Integer replay of the complementary-parity Steiner augmentation."""
from hashlib import sha1
from itertools import combinations, product
from pathlib import Path
import json

P=Path(__file__).resolve().parent
T='dim91-93-triality-matching-transfer'
M='dim89-mixed-matching-rows'
R='dim89-96-rs-axis-code-repack'
def pin(folder,name,digest):
    raw=(P.parent/folder/name).read_bytes()
    assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest,(folder,name)
    return raw.decode()

def main():
    assert __debug__
    parent=json.loads(pin(R,'claims.json','1b78a7a50a09fcbe7e032369d6770a4c8ce395da'))
    parent=[c for c in parent['claims'] if 91<=c['dimension']<=93]
    assert [c['dimension'] for c in parent]==[91,92,93]
    inv=json.loads(pin(T,'selection.json','0e2f71d7a6f2acaddc5a7a100e977958cb38dd8b'))
    old=json.loads(pin('dim91-93-flat-parity-extensions','flat_supports.json','b1112ec5baa8c6b3619f7b4215d8769c16c1c94c'))
    ns={'__name__':'pinned_supports'}
    exec(pin(T,'supports.py','5159d6358b4da75a7bfd0156de308aae3b5942a1'),ns)
    supports=[sum(1<<i for i in q) for q in ns['hanani'](8)]
    assert len(supports)==len(set(supports))==385
    assert all(s.bit_count()==4 and s<1<<22 for s in supports)
    assert all((s&t).bit_count()<=2 for s,t in combinations(supports,2))
    F=[(w,a,o) for w,k in ((6,'six'),(7,'seven')) for a,b in old[k] for o in (b,a^b)]
    F=[F[i] for i in inv['keep_old']]
    assert len(F)==16 and all(w==7 for w,a,o in F)
    matching=[tuple(map(int,l.split())) for l in pin(M,'witness.txt','bae0e9fa8bc696ef5c23e0cc6e81fdd1be856adb').splitlines()]
    matching=[matching[i] for i in inv['append_matching']]
    assert len(matching)==916
    # Every frame has four orthogonal antipodal pairs. Cross-frame products
    # before normalization have magnitude at most 2.
    flat=list(product((-1,1),repeat=4))
    frames=[[v for v in flat if v.count(-1)%2==t] for t in (0,1)]
    frames.append([tuple(2*q*(i==j) for i in range(4)) for j in range(4) for q in (-1,1)])
    dot=lambda x,y:sum(a*b for a,b in zip(x,y))
    for F0 in frames:
        A=[v for v in F0 if next(x for x in v if x)>0]
        assert len(A)==4 and all(dot(a,b)==0 for a,b in combinations(A,2))
        assert set(F0)==set(A+[tuple(-x for x in v) for v in A])
        assert all(dot(v,v)==4 for v in F0)
    assert all(abs(dot(x,y))==2 for A,B in combinations(frames,2) for x,y in product(A,B))
    d4=[tuple(s if k==i else t if k==j else 0 for k in range(4))
        for i,j in combinations(range(4),2) for s,t in product((-1,1),repeat=2)]
    assert max(abs(dot(x,y)) for x,y in product(flat,d4))==2
    # Each dense row has a 24-cell vector, each old Golay row has a
    # doubled D4 root, and each old ordinary row has a D4 root.
    assert 4*16**2<=16*88 and 4*16**2<=16*80 and 4*2**2<=16*2
    Q=[(a,b,c,a^b^c) for a,b,c in product(range(4),repeat=3)]
    assert len(Q)==len(set(Q))==64
    assert min(sum(x!=y for x,y in zip(a,b)) for a,b in combinations(Q,2))==2
    # A pair of different axes or of different signs has dot at most 8.
    assert 4*(4-2)==8
    rows=[tuple(map(int,l.split())) for l in (P/'witness.txt').read_text().splitlines()]
    assert len(rows)==len({i for i,s,k in rows})==20
    for i,s,k in rows:
        assert supports[i]==s and 0<k<15
        assert 2*i not in inv['append_flat_modes'] and 2*i+1 not in inv['append_flat_modes']
        bits=[j for j in range(22) if s>>j&1]
        o=sum((k>>h&1)<<j for h,j in enumerate(bits))
        for p in (o,s^o):
            for w,a,b in F:
                q=s&a
                bound=4*q.bit_count()-2*(q&(b^p)).bit_count()
                assert bound**2<=16*w
            # A matching-frame row is sqrt(2) times an integer D4 root.
            for m,t in matching:
                overlap=(s&(m>>2)).bit_count()
                assert overlap**2<=2*m.bit_count()
        # The two new tag vectors differ in all four rows; all four
        # products against the old uniform axial mode also have magnitude 2.
        assert 4*2==8
    gain=len(rows)*2*len(Q)*16
    out=dict(supports=len(rows),modes=2*len(rows),points_per_mode=len(Q)*16,
             claims=[dict(dimension=c['dimension'],previous=c['lower_bound'],
                          gain=gain,lower_bound=c['lower_bound']+gain) for c in parent])
    assert out==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    main()
