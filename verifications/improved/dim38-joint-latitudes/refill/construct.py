"""Reconstruct the Fano-grid and cube-edge certificate without an optimizer."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import sys

HERE = Path(__file__).resolve().parent

def build():
    recipe = json.loads((HERE/'recipe.json').read_text())
    raw = (HERE.parent/'source/certificate.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest() == recipe['source_certificate_sha256']
    c = json.loads(raw)
    spec = importlib.util.spec_from_file_location('historical38', HERE.parent/'source/verify.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    A, _ = old.leech()
    lookup = {tuple(map(int, row)): i for i, row in enumerate(A)}
    left, right = recipe['left_fano_lines'], recipe['right_fano_lines']
    for offset, lines in ((0, left), (7, right)):
        assert len(lines) == len({tuple(x) for x in lines}) == 7
        assert all(len(x) == len(set(x)) == 3 and all(offset <= j < offset+7 for j in x) for x in lines)
        assert all(len(set(x) & set(y)) == 1 for x,y in itertools.combinations(lines,2))
    pairs = recipe['body_pairs']
    assert len(pairs) == 7 and sorted(sum(pairs, [])) == list(range(14))
    assert recipe['cube_body_coordinates'] == [14,15,16]
    assert recipe['cube_special_positions'] == [[0,1,2],[7,8,9]]
    c['format'] = 'fano_grid_cube_latitudes_v1'
    c['doubled_head_lines'] = []
    c['body_line_indices'] = []
    for i in range(7):
        for j in range(7):
            support = left[i] + right[j]
            a,b = pairs[(i+j)%7]
            for mask in range(32):
                z = [0]*14
                z[support[0]] = 2
                parity = 1
                for k in range(1,6):
                    sign = 1-2*((mask>>(k-1))&1)
                    z[support[k]] = 2*sign
                    parity *= sign
                u = [0]*24
                u[a],u[b] = 4,4*parity
                c['doubled_head_lines'].append(z)
                c['body_line_indices'].append(lookup[tuple(u)])
    c['additional_latitude_points'] = []
    for half in range(2):
        special = recipe['cube_special_positions'][half]
        other = [t for t in range(half*7,half*7+7) if t not in special]
        for i in range(3):
            j,k = (i+1)%3,(i+2)%3
            support = [t for t in range(half*7,half*7+7) if t != special[i]]
            for mask in range(64):
                y = [0]*14
                for t,s in enumerate(support):
                    y[s] = 1-2*((mask>>t)&1)
                parity = 1
                for t in other:
                    parity *= y[t]
                u = [0]*24
                u[14+j],u[14+k] = 4*parity*y[special[j]],-4*parity*y[special[k]]
                sign = 1 if next(v for v in u if v) > 0 else -1
                c['additional_latitude_points'].append(dict(
                    body_line_index=lookup[tuple(sign*v for v in u)],
                    body_sign=sign,doubled_head=[2*v for v in y]))
    c['removed_cap_body_lines'] = sorted(set(c['body_line_indices']) |
        {e['body_line_index'] for e in c['additional_latitude_points']})
    assert len(c['removed_cap_body_lines']) == 20
    c['claimed_size'] = 589680-6*20+13710+1568+28+4*1568+384
    c['previous_local_size'] = 610394
    c['coordinate_rule'] = 'Grid families (+/-u,+/-2Z)/sqrt(128); individual points (body_sign*u,2Z)/sqrt(128); delete all six caps of each listed body line.'
    return (json.dumps(c,indent=2)+'\n').encode()

if __name__ == '__main__':
    output = build()
    if sys.argv[1:] == ['--check']:
        assert output == (HERE/'certificate.json').read_bytes(), 'Deterministic reconstruction'
        print('PASS: byte-for-byte reconstruction, 611522 points.')
    else:
        assert sys.argv[1:] == [], 'Usage: construct.py [--check]'
        (HERE/'certificate.json').write_bytes(output)
