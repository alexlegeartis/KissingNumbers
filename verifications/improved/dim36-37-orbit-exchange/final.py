#!/usr/bin/env python3
"""Derive the compact exchange counts; run verify.py for the complete proof."""
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'dim36-37-coxeter-todd'


def results():
    if not __debug__:raise RuntimeError('Run without -O')
    c=json.loads((HERE/'certificate.json').read_text())
    for name,digest in c['source_sha256'].items():
        assert sha256((SOURCE/name).read_bytes()).hexdigest()==digest
    original=json.loads((SOURCE/'completion/certificate.json').read_text())
    preceding=json.loads((SOURCE/'triangle_switch/proposal.json').read_text())
    old=set(map(tuple,original['phase_orbits']))
    old.difference_update(map(tuple,preceding['remove_phase_orbits']))
    old.update(map(tuple,preceding['add_phase_orbits']))
    removed=set(map(tuple,c['remove_phase_orbits']));added=set(map(tuple,c['add_phase_orbits']))
    assert len(removed)==len(c['remove_phase_orbits']) and removed<=old
    assert len(added)==len(c['add_phase_orbits']) and not(added&old)
    base=2268+672*324+756*243
    previous36=base+3*len(old)
    previous37=previous36+6*len(original['lift_phase_orbits'])
    gain=3*(len(added)-len(removed))
    assert gain==c['gain_each']==9
    assert (previous36,previous37)==(c['previous_points_36'],c['previous_points_37'])
    assert (previous36+gain,previous37+gain)==(c['points_36'],c['points_37'])
    return [{'dimension':36,'total':previous36+gain,'previous':previous36,'gain':gain},
            {'dimension':37,'total':previous37+gain,'previous':previous37,'gain':gain}]


if __name__=='__main__':print(json.dumps(results(),indent=2))
