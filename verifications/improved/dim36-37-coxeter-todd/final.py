#!/usr/bin/env python3
"""Count the current triangle-switch configuration; run verify.py for the proof."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def results():
    source=json.loads((ROOT/'completion/certificate.json').read_text())
    meta=json.loads((ROOT/'triangle_switch/certificate.json').read_text())
    raw=(ROOT/'triangle_switch/proposal.json').read_bytes()
    proposal=json.loads(raw)
    assert meta['format']=='ct_triangle_switch_record_v1'
    assert hashlib.sha256(raw).hexdigest()==meta['proposal_sha256']
    assert hashlib.sha256((ROOT/'completion/certificate.json').read_bytes()).hexdigest()==proposal['source_sha256']==meta['source_sha256']
    original=set(map(tuple,source['phase_orbits']))
    removed=set(map(tuple,proposal['remove_phase_orbits']))
    added=set(map(tuple,proposal['add_phase_orbits']))
    assert len(removed)==len(proposal['remove_phase_orbits']) and removed<=original
    assert len(added)==len(proposal['add_phase_orbits']) and not(added&original)
    flips=proposal['flipped_triangles']
    assert len(flips)==len(set(flips))==meta['flipped_triangles']
    assert all(type(t) is int and 0<=t<61236 for t in flips)
    base=2268+672*324+756*243
    n36=base+3*(len(original)-len(removed)+len(added))
    n37=n36+6*len(source['lift_phase_orbits'])
    assert source['lift_heads']==meta['lift_heads']==[2,-2]
    assert (n36,n37)==(meta['points_36'],meta['points_37'])
    return [{'dimension':36,'total':n36},{'dimension':37,'total':n37}]

if __name__=='__main__':
    for row in results():print(json.dumps(row))
