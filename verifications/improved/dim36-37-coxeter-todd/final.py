#!/usr/bin/env python3
"""Count the phase orbits and latitude layers in the current certificate.

Run verify.py to prove the shell minima and every distance inequality.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def results():
    cert = json.loads((ROOT / 'completion/certificate.json').read_text(encoding='utf-8'))
    assert cert['format'] == 'k12_phase_and_lift_orbits_v2'
    assert cert['squared_norm'] == 12 and cert['lift_heads'] == [2, -2]
    # The fixed base has 2,268 + 672*324 + 756*243 points. The verifier
    # independently enumerates its cosets and derives every multiplicity.
    total36 = 2268 + 672*324 + 756*243 + 3*len(cert['phase_orbits'])
    total37 = total36 + 2*3*len(cert['lift_phase_orbits'])
    assert (total36, total37) == (cert['points_36'], cert['points_37'])
    return [{'dimension': 36, 'total': total36},
            {'dimension': 37, 'total': total37}]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
