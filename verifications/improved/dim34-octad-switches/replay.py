"""Reconstruct the 34D witness and dispatch the established exact replays."""
from copy import deepcopy
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path
from tempfile import TemporaryDirectory

P = Path(__file__).resolve().parent
PRIOR = P.parent / 'dim32-35-signed-layers' / 'rotated34'
DIM = 34


def support(entries):
    assert len(entries) == len(set(entries)) == 8 and all(type(i) is int and 0 <= i < DIM for i in entries)
    return sum(1 << i for i in entries)


def verify(expanded=False, intermediate=False):
    if not __debug__:
        raise RuntimeError("Run without -O")
    raw = (PRIOR / 'certificate.json').read_bytes()
    assert sha256(raw).hexdigest() == "36a9017be6de898d16bdfdad1b686db1d991a9be8a77a804134d770794a5f7ea"
    old = json.loads(raw)
    c = deepcopy(old)
    witness = json.loads((P / 'witness.json').read_text())
    assert len(witness['octads']) == 4 and len(witness['stars']) == 7
    w = c['full_even_supports']
    seen = set()
    for i, before, after in witness['octads'][:2 if intermediate else 4]:
        assert type(i) is int and 0 <= i < len(w) and i not in seen and w[i] == support(before)
        w[i] = support(after)
        seen.add(i)
    assert len(set(w) - set(old['full_even_supports'])) == len(seen)
    stars = []
    for (center, *leaves), parity in witness['stars'][:6 if intermediate else 7]:
        assert len([center, *leaves]) == len(set([center, *leaves])) == 5
        assert all(type(i) is int and 0 <= i < DIM for i in [center, *leaves])
        assert type(parity) is int and parity in (0, 1)
        for a, signs in product((-4, 4), product((-2, 2), repeat=4)):
            if signs.count(-2) % 2 != parity:
                continue
            row = [0]*DIM
            row[center] = a
            for i, s in zip(leaves, signs):
                row[i] = s
            stars.append(row)
    c['stars'] = stars
    assert len(stars) == len(set(map(tuple, stars))) == (96 if intermediate else 112)
    heads = c['heavy_numerators_over_85']
    removed = set()
    for i,j in combinations(range(DIM), 2):
        for a,b in product((-1,1), repeat=2):
            if any(a*x[i]+b*x[j] > 4 for x in stars) or any(4*(a*h[i]+b*h[j]) > 1700 for h in heads):
                removed.add((i,j,a,b))
    assert all((i,j,-a,-b) in removed for i,j,a,b in removed)
    lines = sorted({(i,j,a*b) for i,j,a,b in removed})
    assert len(removed) == 2*len(lines)
    c['removed_root_lines'] = lines
    c['counts'].update(stars=len(stars), roots=2*DIM*(DIM-1)-len(removed))
    c['points'] = sum(c['counts'].values())
    expected = 391000 if intermediate else 391004
    assert c['points'] == expected
    if not intermediate:
        before = set(map(tuple, old['stars']))
        after = set(map(tuple, stars))
        assert (len(after-before), len(before-after)) == (64,16)
        prior_lines = set(map(tuple, old['removed_root_lines']))
        new_lines = set(lines)
        assert (2*len(new_lines-prior_lines),2*len(prior_lines-new_lines)) == (48,12)
        assert len(lines) == 46 and c['counts']['roots'] == 2152
    c['comparison'] = {'preceding_points':390992, 'gain':expected-390992}
    path = PRIOR / ('verify_coordinates.py' if expanded else 'verify.py')
    spec = importlib.util.spec_from_file_location("exact34_layer_proof", path)
    proof = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(proof)
    with TemporaryDirectory() as work:
        p = Path(work) / 'certificate.json'
        p.write_text(json.dumps(c, sort_keys=True) + '\n')
        options = dict(expected_points=expected, expected_lines=len(lines))
        if not expanded:
            options['expected_stars'] = len(stars)
        out = proof.verify(p, **options)
    assert out['points'] == expected and out['all_checks_pass']
    out.pop('elapsed_seconds', None)
    return out
