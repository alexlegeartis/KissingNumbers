"""Derive the cardinality from the witness; verify.py proves its geometry."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def results():
    cert = json.loads((ROOT / 'source/certificate25.json').read_text())
    circles = sum(len(cl['owners']) for cl in cert['classes'])
    removed = {tuple(u) for cl in cert['classes'] for u in cl['owners']}
    removed.update(tuple(a + b for a, b in zip(u, cl['shift']))
                   for cl in cert['classes'] for u in cl['owners'])
    removed.update(map(tuple, cert['local_owners'] + cert.get('new_owners', [])))
    small = len(cert['small'])
    total = 196560 - len(removed) + 3 * circles + small
    if (circles, len(removed), small, total) != (971, 1991, 103, 197585):
        raise ValueError('Unexpected witness counts')
    if total != cert['points']:
        raise ValueError('Certificate point count disagrees with its data')
    return [{'dimension': 25, 'total': total}]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
