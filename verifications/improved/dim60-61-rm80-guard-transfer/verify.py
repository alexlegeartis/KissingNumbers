#!/usr/bin/env python3
"""Transfer exact signed support banks into the published 60/61 architectures."""
from collections import defaultdict,Counter
from contextlib import redirect_stdout
from io import StringIO
from itertools import combinations
from pathlib import Path
import runpy,sys,json,gzip,base64
from hashlib import sha256
assert __debug__
here=Path(__file__).resolve().parent
src=here.parent/'dim60-61-rm80-puncture'
PINS={'verify_extra.py': 'a90c346b862c6f2c4cfcefa1e12c13ba58e4443145c12d2c802c28e290c63b41', 'extra61.json.gz.b64': 'a99b0fe9f4f91ff0407e7a66c820ab0f2243477d1c84d66159aedb1f308d42e7', 'verify.py': 'ebee4932f44f7f3e15eadeadcc9f7bf4c4eea206aa284621e4cd3a4bb5e79d67', 'extra60.json.gz.b64': '8420af6106dcc6a6a2bda865593813e14e36f76dc547bac2accfcce320d05ac2'}
assert all(sha256((src/p).read_bytes()).hexdigest()==v for p,v in PINS.items())
with redirect_stdout(StringIO()):
    V=runpy.run_path(str(src/'verify_extra.py'))
marks=V['marks']
R=V['rows']
bounds={}
for n,base in [(60,68033576),(61,118639402)]:
    modes={0:0} if n==60 else {0:0,1:1,1<<20:-1}
    D=json.loads(gzip.decompress(base64.b64decode((src/f'extra{n}.json.gz.b64').read_bytes())))
    old={s for s in V['oldAll'] if (s&marks) in modes}|set(D['extra'])
    L=defaultdict(list)
    for a,b,g,t in D['guards']:
        L[a].append((g,t^1));L[b].append((g,t))
    B=gzip.decompress(base64.b64decode((here/f'transfer{n}.txt.gz.b64').read_bytes())).splitlines()
    p,q=map(int,B[0].split())
    a=list(map(int,B[1:1+p]));assert len(a)==len(set(a))==p and not set(a)&old
    N=set(a);S=old|N
    assert len(B)==1+p+q
    guards={}
    for row in B[1+p:]:
        x,y,g,t=map(int,row.split())
        assert x in N and y in S and x!=y and t in (0,1)
        assert g.bit_count()==8 and g&(x&y)==g
        k=(min(x,y),max(x,y));assert k not in guards
        guards[k]=g
        L[x].append((g,t^1));L[y].append((g,t))
    for v in S:
        assert v.bit_count()==16 and v&marks in modes
        assert all(not ((v&w).bit_count()&1) for w in R)
        H={0:0}
        for g,t in [(v,int(not v&marks))]+L[v]:
            if g in H: assert H[g]==t
            else: H.update({h^g:r^t for h,r in list(H.items())})
        assert len(H) in (2,4,8,16,32) and all(g in (0,v) or g.bit_count()==8 for g in H)
    def geo(x,y):
        t=(x&y&~marks).bit_count()+modes[x&marks]*modes[y&marks]
        if t<=8:return 0
        assert t<=10 and (min(x,y),max(x,y)) in guards
        return 1
    actual=sum(geo(x,y) for x in N for y in old)+sum(geo(x,y) for x,y in combinations(a,2))
    assert actual==q
    cnt=Counter(x&marks for x in N)
    gain=sum(v*(1024 if t else 2048) for t,v in cnt.items())
    assert gain==({60:688128,61:210944}[n])
    bounds[str(n)]=base+gain
    print(json.dumps(dict(dimension=n,added_banks=dict(cnt),guards=q,gain=gain,lower_bound=base+gain)))
assert {'status':'PASS','lower_bounds':bounds}==json.loads((here/'claims.json').read_text())
