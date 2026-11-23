#!/usr/bin/env python3
"""Exact spike-amplified DG layers for the transferred RM80 60/61 stocks."""
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations_with_replacement
from io import StringIO
from contextlib import redirect_stdout
from hashlib import sha256
import json,runpy,sys
from dg import code
assert __debug__
H=Path(__file__).resolve().parent
parent=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else H.parent
source=next(p for p in (parent/'dim60-61-rm80-guard-transfer',parent/'release601') if p.is_dir())
pins={'claims.json':'5a423875aba574e710c78d40630825c3d53e41daebacf6f138e9ce0a60a55b5f',
      'verify.py':'b21509276df4a6cb5d52e0bf00efddac066fc53665b2c1c47598d7db9f2330a8',
      'transfer60.txt.gz.b64':'7fb3b5ed5b309bc87f05df08985e45e0f9d14c240835f4404b6912b748dcf65e',
      'transfer61.txt.gz.b64':'8d1cefaa459b09f165eb161ee123739b126dde761792e5aeb121c8e833cf2a08'}
assert all(sha256((source/f).read_bytes()).hexdigest()==v for f,v in pins.items())
with redirect_stdout(StringIO()):
    sys.argv=[str(source/'verify.py')]
    B=runpy.run_path(str(source/'verify.py'))
assert B['bounds']=={'60':68721704,'61':118850346}
marks=(0,20,41,61);mask=sum(1<<i for i in marks)
assert B['marks']==mask
V=B['V']['v'];assert V['x']['poles']==120 and V['y']['poles']==122
C,A=code()
assert len(C)==1024 and len(A)==128
cube=V['C'];Q=V['R']
assert len(cube)==80 and len(Q)==22
assert len(set(cube[:60]))==60
counts=Counter(v&mask for v in A);assert len(counts)==8 and set(counts.values())=={16}
assert all(v.bit_count()%2==0 for v in counts)
def pat(t):return sum(1<<v for i,v in enumerate(marks)if t>>i&1)
def gate(x,u,v):return x<=0 or 4*x*x<=u*v
def key(t,h):return (t,F(h))
config={
60: ([(2,F(11,2),F(0)),(3,F(11,2),F(0))],[(0,F(0)),(1,F(0))],[((),0,1)]),
61: ([(7,F(29,5),-F(7,4)),(5,F(29,5),F(7,4)),(6,F(43,8),-F(3,5))],
     [(0,F(1,2)),(1,F(2)),(2,-F(1,2)),(3,-F(2))],
     [((),0,1),((0,),1,0),((1,),-1,0)]),
}
assert Counter((q&mask).bit_count()%2 for q in C)=={0:512,1:512}
claims={}
for n,(donors,base,sparse) in config.items():
    assert len({t for t,a,u in donors})==len(donors) and not {s for s,h in base}&{t for t,a,u in donors}
    assert {t.bit_count()%2 for t,a,u in donors}=={0,1}
    for t,a,u in donors:
        N=59+a*a+u*u
        assert a>1
        for s,h in base:
            # These labels are absent from the old dense carrier, so the
            # original full words are distinct even for one cubic coset.
            d=16-(s^t).bit_count()
            assert gate(60-2*d+a-1+u*h,N,60+h*h)
        for mk,h,parity in sparse:
            agree_parity=sum(t>>i&1 for i in mk)%2
            # Full support orthogonality to the mother code is inherited,
            # so a parity disagreement forces at least one physical
            # mismatch, costing at least two even in the spiked body site.
            dot=16-len(mk)+a-1+u*h-2*(agree_parity!=parity)
            assert gate(dot,N,16)
        assert gate(3+a,N,4)
        if n==61:
            assert gate(2+a+abs(u),N,4) and gate(abs(u),N,1)
        assert not gate(a,N,1)  # remove all 120 ordinary poles
    for (s,a,h),(t,b,u) in combinations_with_replacement(donors,2):
        N=59+a*a+h*h;M=59+b*b+u*u
        d=(s^t).bit_count()
        # Same spike: use the DG d24 subcode. Distinct spike positions use
        # distinct mother cubic cosets and the proven parent d16.
        assert gate(60-2*(24-d)+a*b-1+h*u,N,M)
        assert gate(60-2*(16-d)+a+b-2+h*u,N,M)
    layer=[sum(counts.get((q^c^pat(t))&mask,0)for c in cube[:60]for q in C)for t,a,h in donors]
    assert all(x==491520 for x in layer)
    add=sum(layer)
    count=B['bounds'][str(n)] + add -120
    assert add=={60:983040,61:1474560}[n]
    claims[str(n)]=count
    print(json.dumps(dict(status='PASS',dimension=n,prior=B['bounds'][str(n)],added_spikes=add,
               dropped_ordinary_poles=120,net_gain=add-120,lower_bound=count,
               labels=[t for t,a,h in donors],norms_squared=[str(59+a*a+h*h)for _,a,h in donors])))
assert {'status':'PASS','lower_bounds':claims}==json.loads((H/'claims.json').read_text())
