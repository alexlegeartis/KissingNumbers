#!/usr/bin/env python3
"""Ceilings for families restricted to the shipped P48 base classes.

For a fixed base of C lines, any disjoint family of its automorphism images
has S(lambda) <= lambda*C. This does not bound larger base classes or other
constructions. The 7,387-line October class supersedes the old 7,069-line
base; neither is a proved maximum independent set of the mother graph.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from october_family import families

ERS0 = {48: 24, 52: 21, 56: 23, 60: 25, 61: 25, 62: 26, 63: 27}


def main():
    rows = json.loads((ROOT / 'rows49-63.json').read_text())
    bases = families()
    maximum = max(int(sizes.max()) for sizes in bases.values())
    assert maximum == 7387
    print(__doc__)
    print('dim  current cap total  fixed-base ceiling  ERS level 0')
    for row in rows:
        dim, count = row['dim'], row['lam']
        poles = row['gain'] - 2*row['S']
        ceiling = 52416000 + 2*count*maximum + poles
        assert row['value'] <= ceiling
        floor = max([1 << k for n,k in ERS0.items() if n <= dim] or [0])
        print(dim, row['value'], ceiling, floor)
    print('PASS: fixed-base ceilings checked. Larger classes remain open.')


if __name__ == '__main__':
    main()
