#!/usr/bin/env python3
"""Exact simultaneous support/phase audit; the solver is not used."""
from collections import defaultdict
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from itertools import combinations
from fractions import Fraction
import base64,gzip,json,runpy
assert __debug__
root=Path(__file__).resolve().parent
with redirect_stdout(StringIO()):
    v=runpy.run_path(str(root/'verify.py'))
M=(0,20,41,61)
marks=sum(1<<p for p in M)
oldAll=v['V']
rows=v['R']+v['C']
assert len(oldAll)==37996 and len(set(oldAll))==len(oldAll)
for n in (60,61):
    modes={0:0} if n==60 else {0:0,1:1,1<<20:-1}
    prior=v['x' if n==60 else 'y']['total']
    old={s for s in oldAll if (s&marks) in modes}
    D=json.loads(gzip.decompress(base64.b64decode((root/f"extra{n}.json.gz.b64").read_bytes())))
    X=D['extra']; N=set(X)
    assert len(X)==len(N) and not old&N
    all_=old|N
    assert all(s.bit_count()==16 and s&marks in modes for s in all_)
    assert all((s&g).bit_count()%2==0 for s in all_ for g in rows)
    syn=defaultdict(list)
    guard={}
    # Apply the inherited odd ordinary banks and even forced marked banks.
    # The marked parities are even because both unique facing labels have
    # a negative sign bit on their removed mark.
    for s in all_: syn[s].append((s,int(not(s&marks))))
    for a,b,g,t in D['guards']:
        assert a in N and b in all_ and a!=b and t in (0,1)
        assert g.bit_count()==8 and g&(a&b)==g
        k=tuple(sorted((a,b)));assert k not in guard
        guard[k]=g
        syn[a].append((g,1-t)); syn[b].append((g,t))
    # Treat a deleted mark as a ZERO SIGN dummy coordinate, replacing
    # it physically by the fixed directed height. Every nontrivial full
    # check row other than total parity has weight exactly eight. Hence
    # the affine constant-weight check space extends to the 5-check
    # extended Hamming dual by repeatedly bisecting equal-size classes.
    for s,L in syn.items():
        A={0:0}
        for g,p in L:
            if g in A:
                assert A[g]==p
            else:
                A.update({a^g:b^p for a,b in list(A.items())})
        assert len(A) in (2,4,8,16,32)
        assert all(g in (0,s) or g.bit_count()==8 for g in A)
        # On a marked support, shortening at the one dummy costs exactly
        # half the signs: a singleton parity is not in the completed dual.
        assert not s&marks or (s&marks).bit_count()==1
    def test(a,b):
        cap=((a&b)&~marks).bit_count()+modes[a&marks]*modes[b&marks]
        if cap<=8:return 0
        assert cap<=10
        s=guard[tuple(sorted((a,b)))]
        assert ((s&~marks)&(a&b))==(s&~marks) and (s&~marks).bit_count()>=7
        return 1
    hits=sum(test(a,b) for a in N for b in old)+sum(test(a,b) for a,b in combinations(X,2))
    # A new active guard contributes at least one ordinary sign mismatch;
    # dummy fixed-zero sign cannot absorb it.
    sizes={t:sum((s&marks)==t for s in X) for t in modes}
    stock=sum(z* (1024 if t else 2048) for t,z in sizes.items())
    # Dense, planes and poles use the same per-mode comparisons as base.
    assert prior+stock==({60:68033576,61:118639402}[n])
    print(json.dumps(dict(dimension=n,added_supports=sizes,total_guards=len(guard),
                          dangerous_pairs=hits,added_points=stock,lower_bound=prior+stock)))
