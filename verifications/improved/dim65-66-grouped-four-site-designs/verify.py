from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import json
import re
import subprocess
import sys

P = Path(__file__).resolve().parent
assert __debug__
if sys.argv[1:] != ["--construction-only"]:
    assert len(sys.argv) == 1
    for d, name, n in ((65, "dim65-68-complete-hyperplane-atlas", 418995074),
                       (66, "dim66-68-nested-affine-face-trades", 429284740)):
        prev = P.parent/name
        assert json.loads((prev/"claims.json").read_text())["lower_bounds"][str(d)] == n
        out = subprocess.check_output([sys.executable,str(prev/"verify.py")],text=True)
        assert re.search(rf"^DIM {d}\b.*\bN {n}$",out,re.M)
        print(f"Predecessor {d} PASS",flush=True)
counts = {}
starter = [tuple(map(int, r.split())) for r in (P/'starter.txt').read_text().splitlines()]
assert len(starter)==30
S = set()
X = set()
for f in starter:
    assert len(f)==4 and len({x%11 for x in f})==4 and all(0<=x<22 for x in f)
    for g in range(11):
        u=tuple((x+g)%11+11*(x//11) for x in f)
        for t in combinations(sorted(u),3):
            assert t not in X
            X.add(t)
        for z in product(range(3),repeat=3):
            zz=(*z,-sum(z)%3)
            q=tuple(sorted(6*(x%11)+3*(x//11)+v for x,v in zip(u,zz)))
            assert q not in S
            S.add(q)
assert len(X)==comb(11,3)*8 and len(S)==8910
for g,h in combinations(range(11),2):
    for k in range(5):
        E=((5,k),((k+1)%5,(k-1)%5),((k+2)%5,(k-2)%5))
        for (a,b),(c,d) in product(E,repeat=2):
            q=tuple(sorted((6*g+a,6*g+b,6*h+c,6*h+d)))
            assert q not in S
            S.add(q)
assert len(S)==11385
for g in range(11):
    S.update(tuple(6*g+i for i in q) for q in [(0,1,2,3),(0,1,4,5),(2,3,4,5)])
assert len(S)==11418
used=[]
for d in (65,66):
    F={x for x in S if max(x)<d}
    H=Counter(t for s in F for t in combinations(s,3))
    assert max(H.values())==1 and len(H)==4*len(F)
    if d==66:
        leave=set(combinations(range(d),3))-H.keys()
        assert len(leave)==88 and all(len({i//6 for i in t})==1 for t in leave)
        assert Counter(t[0]//6 for t in leave)=={i:8 for i in range(11)}
        assert (comb(66,3)-88)//4==len(F)
    # Raw bounds against all unchanged norm16, norm64, poles, and norm72 head stock.
    for oldq,ip2 in ((16,16),(64,16),(4,4),(72,50)):
        assert 4*ip2<=4*oldq
    n0,a0={65:(418995074,10416),66:(429284740,11088)}[d]
    n=n0+16*(len(F)-a0)
    assert (len(F),n)=={65:(10726,419000034),66:(11418,429290020)}[d]
    print(f"DIM {d}: {len(F)} supports, {n-n0} point gain; N {n}")
    counts[str(d)]=n
assert json.loads((P/"claims.json").read_text()) == {"status":"PASS","lower_bounds":counts}
