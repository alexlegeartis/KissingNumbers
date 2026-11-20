"""Rebuild the explicit 33D certificate and check it in exact arithmetic."""
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path
from tempfile import TemporaryDirectory

P = Path(__file__).resolve().parent
Q = P.parent / 'dim32-35-signed-layers' / 'nonantipodal33'


def mask(indices):
    assert len(indices) == len(set(indices)) == 8
    assert all(type(j) is int and 0 <= j < 33 for j in indices)
    return sum(1 << j for j in indices)


def replay(expanded=False):
    if not __debug__:
        raise RuntimeError("Run without -O")
    raw = (Q / 'certificate33_370408_explicit.json').read_bytes()
    assert sha256(raw).hexdigest() == "5686dfcd2b891bf6d9be6b01013d8205774aa1a19c02a2c45bdfd86194314a60"
    c = json.loads(raw)
    w = json.loads((P / 'witness.json').read_text())
    assert c['points'] == 370408 and len(w['octads']) == len(w['stars']) == 2
    changed = set()
    for i, old, new in w['octads']:
        assert type(i) is int and 0 <= i < 1837 and i not in changed
        assert c['full_even_supports'][i] == mask(old) != mask(new)
        c['full_even_supports'][i] = mask(new)
        changed.add(i)
    stars = c['star_signed_points']
    assert len(stars) == 180
    for i,*L in w['stars']:
        assert len(L) == 4 and len(set([i,*L])) == 5
        assert all(type(j) is int and 0 <= j < 33 for j in [i,*L])
        for a, leaves in product((-4,4),product((-2,2),repeat=4)):
            if leaves.count(-2) % 2:
                continue
            x = [0]*33
            x[i] = a
            for j,b in zip(L,leaves):
                x[j] = b
            stars.append(x)
    assert len(stars) == 212 and len(set(map(tuple, stars))) == 212
    old_removed = set(map(tuple,c['removed_root_signed_points']))
    removed = []
    for i,j in combinations(range(33),2):
        for a,b in product((-4,4),repeat=2):
            if any(a*x[i]+b*x[j]>16 for x in stars):
                v = [0]*33
                v[i],v[j] = a,b
                removed.append(v)
    assert len(old_removed) == 156 and old_removed <= set(map(tuple,removed))
    assert len(removed) == 180
    c['removed_root_signed_points'] = removed
    c['points'] = 370416
    c['construction'] = {'preceding_points':370408, 'gain':8, 'octad_switches':2}
    spec = importlib.util.spec_from_file_location(
        "exact33", Q / ('verify_coordinates.py' if expanded else 'verify.py'))
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    with TemporaryDirectory() as directory:
        path = Path(directory) / "certificate.json"
        path.write_text(json.dumps(c, sort_keys=True) + '\n')
        check = verifier.verify(path)
    assert check['points'] == 370416 and check['all_checks_pass']
    if not expanded:
        assert check['unpaired_points'] == 500
    check.pop('seconds',None)
    return check
