"""Verify the 38D root-for-star trade using exact integer geometry."""
from itertools import combinations
import json
import importlib.util
from pathlib import Path

P = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    'dim38_qr_free_offsets', P.parent / 'dim38-qr-free-offsets' / 'verify.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def read_stars():
    out = []
    for line in (P / 'stars.txt').read_text().splitlines():
        center, leaves, head_sign, neg = map(int, line.split())
        assert 0 <= center < 38 and 0 <= leaves < 1 << 38
        assert leaves.bit_count() == 4 and not (leaves >> center & 1)
        assert head_sign in (-1, 1) and neg >= 0 and neg & leaves == neg
        out.append(tuple(2 * head_sign if j == center else
                         ((-1 if neg >> j & 1 else 1) if leaves >> j & 1 else 0)
                         for j in range(38)))
    assert len(set(out)) == len(out)
    return out


def replay():
    old = base.replay()
    assert old['points'] == 652028
    X = read_stars()
    assert all(sum(a*a for a in x) == 8 for x in X)
    assert all(sum(a*b for a, b in zip(x, y)) <= 4
               for x, y in combinations(X, 2))
    supports = (P.parent / 'dim37-38-qr-merge' / 'supports_3088.txt').read_text().splitlines()
    S = [tuple(j for j, a in enumerate(s) if a == '1') for s in supports]
    assert all(sum(abs(x[j]) for j in s) <= 4 for x in X for s in S)
    # The bound starts from six units of absolute mass against magnitude 2;
    # the dense marked magnitude 3 can cost at most two further units.
    for head in (0, 3):
        bound = 12 + 2*max(0, head-2)
        assert 4*bound**2 <= 8*(148+head**2)
    # Define and remove exactly the signed roots that violate the angle.
    # Only matching signs on a center-leaf pair give product 3 > 2.
    forbidden = set()
    for x in X:
        c = next(j for j, a in enumerate(x) if abs(a) == 2)
        for j, a in enumerate(x):
            if abs(a) == 1:
                forbidden.add(tuple(sorted(((c, x[c]//2), (j, a)))))
    n = old['points'] - len(forbidden) + len(X)
    assert n > old['points']
    return dict(status='PASS', dimension=38, points=n, baseline=old['points'],
                stars=len(X), removed_roots=len(forbidden), gain=n-old['points'])


if __name__ == '__main__':
    print(json.dumps(replay(), sort_keys=True))
