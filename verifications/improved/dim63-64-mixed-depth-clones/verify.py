#!/usr/bin/env python3
"""Exact mixed-amplitude dense-owner clones in 63D and 64D."""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json, importlib.util
HERE=Path(__file__).resolve().parent
def verify():
    assert __debug__
    f=HERE.parent/'dim63-64-dense-owner-clones/verify.py'
    assert sha256(f.read_bytes()).hexdigest()=='232836145af472c07dc464b58d6a2041a7ecbf5f63986739b04736a3a6ecb480'
    p=importlib.util.spec_from_file_location('prior',f)
    old=importlib.util.module_from_spec(p);p.loader.exec_module(old)
    parent=old.verify()['dimensions']
    G,J,R=old.stock()
    Z=old.load('dim63-65-affine-reconstruction/section63.json',
               'cdc01424d1789116ac83c5bb6509d88269bfc16b')['octad']
    assert Z.bit_count()==8
    # All points have squared norm 64. h is the heavy unscaled integer,
    # ordinary unscaled coordinates have magnitude 2, and scale is 8/sqrt(N).
    htype=(11,14); N={h:63*4+h*h for h in htype}
    for h in htype:
        # Original sparse W4, W9 envelope, plus retained axial poles.
        for v in (Q(4*(h+6)),Q(8,3)*(h+16),Q(16)):
            assert v*v<=16*N[h]
        # Original W16 needs either an ordinary mismatch on a support
        # containing the mark, or no mark. Both cases pass exactly.
        for v in (Q(2*(h+26)),Q(64)):
            assert v*v<=16*N[h]
        assert (h+66)**2<=16*N[h]
    # Distinct owners on distinct marks: parent d>=16 and two changed
    # signs at different marks yield this exact upper envelope.
    for h in htype:
        for g in htype:
            assert 4*(136+2*h+2*g)**2<=N[h]*N[g]
    # Same owner gets at most one shallow mark. Distinct mark safe
    # with deep+shallow or deep+deep.
    for h,g in ((11,14),(14,14)):
        assert 4*(248-2*h-2*g)**2<=N[h]*N[g]
    out=[]
    for d,pr in zip((63,64),parent):
        assert pr['dimension']==d
        W=set(); used=0; shallow=0;deep=0;S=[];per=[[] for _ in range(64)]
        for ln in (HERE/f'win{d}.txt').read_text().splitlines():
            w,M,L=map(int,ln.split())
            assert 0<=w<1<<64 and 0<M<1<<64 and 0<=L<1<<64
            assert w not in W and M&L==0 and L.bit_count()<=1 and (M|L).bit_count()>=2
            assert d==64 or ((w&Z).bit_count()==4 and not (M|L)&Z)
            v=w
            for row,j in zip(G,J):
                if v>>j&1:v^=row
            assert v==0
            for mask,h in ((M,14),(L,11)):
                for a in old.bits(mask):
                    x=w^(1<<a)
                    assert all(any((u&x).bit_count()%2!=t for u,t in c) for c in R[a])
                    assert d==64 or (x&Z).bit_count()==4
                    S.append((a,x,h));per[a].append((x,h))
            W.add(w);used|=M|L;shallow|=L;deep|=M
        assert shallow&deep==0
        assert len(S)==len(set((a,x) for a,x,h in S))
        same_pairs=0
        for a,P in enumerate(per):
            for (x,h),(y,g) in combinations(P,2):
                v=x^y; da=v>>a&1
                dot=4*(63-2*(v.bit_count()-da))+h*g*(1-2*da)
                assert dot<=0 or 4*dot*dot<=N[h]*N[g]
                same_pairs+=1
        poles=2*used.bit_count()
        total=pr['construction_base']+len(S)-len(W)-poles
        out.append(dict(dimension=d,previous=pr['lower_bound'],stock=pr['construction_base'],
                        lower_bound=total,gain=total-pr['lower_bound'],inserted=len(S),
                        removed_dense=len(W),removed_poles=poles,
                        shallow_clones=sum(h==11 for _,_,h in S),shallow_marks=shallow.bit_count(),
                        exact_same_mark_pairs=same_pairs))
    return dict(dimensions=out)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
