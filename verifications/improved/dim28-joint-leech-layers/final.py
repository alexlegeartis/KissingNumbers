"""Derive the 28D count from the supplied partition and auxiliary points."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def results():
    classes, auxiliary = [], 0
    with (ROOT / 'source/D28-205232.txt').open(encoding='utf-8') as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == 'dimension':
                dimension = int(fields[1])
            elif fields[0] == 'points':
                declared = int(fields[1])
            elif fields[0] == 'class':
                classes.append((len(fields[1].split(',')), int(fields[2])))
            elif fields[0] in ('extra', 'pole'):
                auxiliary += 1
    assert dimension == 28 and classes[0] == (1, 192352)
    assert len(classes) == 9 and sum(n for _, n in classes[1:]) == 4208
    assert all(tails == 3 for tails, _ in classes[1:]) and auxiliary == 256
    total = sum(tails * size for tails, size in classes) + auxiliary
    assert total == declared == 205232
    return [{'dimension': dimension, 'total': total}]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
