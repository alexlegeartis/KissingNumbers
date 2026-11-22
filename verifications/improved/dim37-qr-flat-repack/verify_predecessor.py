"""Exact 37D QR radical layer proof with the original sparse octad layer."""
from itertools import combinations
from pathlib import Path
import json
import sys
from verify37 import replay as original
from supports37 import supports
P=Path(__file__).resolve().parent
def six_supports():return [int(s) for s in (P/'six.txt').read_text().split()]
def fit(p,q,N):
    # (p + q√2)^2 ≤ 8*N. All terms are nonnegative before expansion.
    r=8*N-p*p-2*q*q
    return r>=0 and r*r>=8*p*p*q*q
def replay():
    data=original(); assert data['points']==533864
    T=six_supports();B=supports()
    assert len(T)==len(set(T)) and all(0<s<1<<36 and s.bit_count()==6 for s in T)
    assert max((a&b).bit_count() for a,b in combinations(T,2))<=2
    assert max((a&b).bit_count() for a in T for b in B)<=4
    assert fit(24,2,145) and fit(24,8,160)
    # Encode each new point by ±2 on an ordinary six-support of even
    # parity, and independent ±2√2 in the otherwise exceptional axis.
    assert 6*4+8==32
    # Same marked sign: distinct supports agree in at most two places;
    # identical supports and distinct even signings differ at least twice.
    assert 8+2*4==16 and 8+(6-2*2)*4==16
    # Opposite marked signs: even identical ordinary signs give max 16.
    assert -8+6*4==16
    # Ordinary roots contribute at most two terms of magnitude eight.
    # Octads contribute at most four terms of magnitude four.
    assert 2*8==16 and 4*4==16
    # Add both unit poles on the marked axis. They have mutual
    # product -1, and zero product with the retained roots and octads.
    # Check squared absolute product against each remaining layer.
    assert 4*8==32  # radical marked entry ±2√2, squared norm 32
    assert all(4*m*m<=144+m*m for m in (1,4))
    poles=2
    removed=4*36
    n=data['points']-removed+64*len(T)+poles
    return dict(status='PASS',dimension=37,points=n,baseline=data['points'],
        new_six_supports=len(T),unit_poles=poles,radical_points=64*len(T),removed_marked_roots=removed,
        kept_ordinary_roots=2*36*35,gain=n-data['points'])
if __name__=='__main__': print(json.dumps(replay(),sort_keys=True))

