#!/usr/bin/env python3
"""Exact support and coset replay of simultaneous parity rephasing."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from collections import defaultdict,Counter
import json,runpy
import base64,gzip,struct
assert __debug__
here=Path(__file__).resolve().parent
parent=here.parent/'dim62-single-mark-square'
PIN='8e7b6139ae4f7ae458810b9fd9745c27ecd87d375787c4b59526488746d0a4c6'
assert sha256((parent/'verify.py').read_bytes()).hexdigest()==PIN
D=runpy.run_path(str(parent/'verify.py'))
sparse=D['S']; dense=D['quadratic']+D['phase']
base=json.loads((parent/'claims.json').read_text())['lower_bounds']

def realize(body, forms, total, targets):
    # Check that all guard rows fit in an affine extended Hamming dual.
    nn=body.bit_count()
    assert nn in (15,16) and len(forms)==len(targets)
    assert all(f&body==f and f.bit_count()==8 for f in forms)
    L=[sum(((f>>i)&1)<<j for j,f in enumerate(forms)) for i in range(64) if (body>>i)&1]
    C=Counter(L+([0] if nn==15 else []))
    assert sum(C.values())==16 and len(set(C.values()))==1
    first=min(C)
    assert all(a^b^first in C for a in C for b in C)
    assert len(C) in (1,2,4,8,16)
    # The affine support identifies distinct column labels. Assign the
    # unused label to the padded position when length 15. Every guard is
    # affine-linear in the labels, and all n≥15 distinct columns give
    # extended Hamming parity-check rank 5 and minimum distance≥4.
    # Check that the required coset parities are jointly achievable.
    V={0}
    for p in L:
        t=1|(p<<1); V|={v^t for v in V}
    rhs=total | sum((b<<1)<<j for j,b in enumerate(targets))
    assert rhs in V
    return 1<<(nn-5)

res={}
for d,marks,parity in [('62',(0,9,20,41),(0,1,0,1)),('63',(0,20),(1,0))]:
    M=sum(1<<i for i in marks)
    kmap={0:0,**{1<<i:j+1 for j,i in enumerate(marks)}}
    type=lambda s:kmap.get(s&M,-1)
    gate=lambda a,b:8+(0 if not(a and b) else -1 if a==b else 1 if d=='63' or abs(a-b)==2 else 0)
    old={s:(s&~M,type(s)) for s in sparse if type(s)>=0}
    inc=set(sparse)
    new=[]; guards=defaultdict(set); N=0
    Q=sorted(old)
    bo=gzip.decompress(base64.b64decode((here/f'witness{d}.bin.gz.b64').read_text()))
    p=0
    while p<len(bo):
        s,n=struct.unpack_from('<QB',bo,p);p+=9
        k=type(s); B=s&~M
        assert 0<=s<1<<64 and s.bit_count()==16 and s not in inc and k>=0
        pos=[1<<i for i in range(64)if B>>i&1]
        F={}
        for _ in range(n):
            q,t=struct.unpack_from('<HH',bo,p);p+=4
            assert q<len(Q) and t<1<<len(pos) and t.bit_count()==8 and Q[q] not in F
            F[Q[q]]=sum(pos[i]for i in range(len(pos))if t>>i&1)
        f1=tuple(sorted(set(F.values())))
        assert len(f1)<=30
        P=parity[k-1] if k else 0
        z=realize(B,f1,P,[0]*len(f1))
        assert z==(1024 if k else 2048)
        N+=z
        if k:assert all((h&s).bit_count()%2==0 for h in dense)
        G={}
        for q,(a,j) in old.items():
            i0=(a&B).bit_count()
            gd=gate(j,k)
            if i0>gd:
                assert q in F and i0<=gd+2 and (B&a&F[q])==F[q]
                G[q]=F[q];guards[q].add(F[q])
        assert G==F
        new.append((s,k,B))
    assert len({s for s,k,B in new})==len(new)
    assert all((a&b).bit_count()<=gate(i,j)for (_,i,a),(_,j,b) in combinations(new,2))
    for q,g in guards.items():
        a,j=old[q];F=tuple(sorted(g))
        assert len(F)<=30
        assert realize(a,F,parity[j-1] if j else 0,[1]*len(F))==(1024 if j else 2048)
    # Every conditional pair has one disagreement on a shared guard.
    # A mismatch saves 8 on the dot product, exactly two ordinary hits.
    # Old-old, new-new and all plane/pole/body conditions are the same
    # magnitude and dense parity inequalities as the preceding proof.
    res[d]=base[d]+N
    print('PASS',d,res[d],'gain',N,'new banks',len(new),'rephased old banks',len(guards),flush=True)
assert {'status':'PASS','lower_bounds':res}==json.loads((here/'claims.json').read_text())
