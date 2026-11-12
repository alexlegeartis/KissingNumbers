"""Replay manifest, complete exact geometry and independent changed pairs."""
from pathlib import Path
import hashlib
import json
import sys
from verify_geometry import verify
from verify_delta import replay

ROOT = Path(__file__).resolve().parent


def check_manifest():
    manifest = json.loads((ROOT/'data/manifest.json').read_text())
    actual = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in ROOT.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.name != 'manifest.json'}
    if actual != manifest:
        raise ValueError('Complete file manifest')


if __name__ == '__main__':
    wanted = set(map(int, sys.argv[1:]))
    if wanted - {28, 29}:
        raise ValueError('Supported dimensions are 28 and 29')
    check_manifest()
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        if wanted and record['dimension'] not in wanted:
            continue
        geometry = verify(ROOT/record['certificate'])
        expected = json.loads((ROOT/record['verification']).read_text())
        if {k: v for k, v in geometry.items() if k != 'seconds'} != {
                k: v for k, v in expected.items() if k != 'seconds'}:
            raise ValueError('Complete geometry report')
        delta = replay(record)
        if delta != record['direct_integer_replay']:
            raise ValueError('Independent changed-pair report')
        print(json.dumps(dict(status='PASS', geometry=geometry, delta=delta)), flush=True)
