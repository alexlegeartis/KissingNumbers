#!/usr/bin/env python3
"""Replay the exact 36D phase code and its latitude extension in 37D.

All generated reports go into a temporary copy, preserving the source package.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from final import ROOT, results


def main():
    if not __debug__:
        raise SystemExit('Run without -O: verification assertions must be enabled.')
    manifest = ROOT / 'data/SHA256SUMS.txt'
    listed = {}
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        path = (ROOT / name).resolve()
        assert ROOT in path.parents and path != manifest and name not in listed
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
        listed[name] = digest
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p != manifest}
    assert set(listed) == actual, 'manifest coverage'
    with tempfile.TemporaryDirectory(prefix='kissing36_37_') as directory:
        target = Path(directory) / 'phase'
        shutil.copytree(ROOT / 'phase', target, ignore=shutil.ignore_patterns('__pycache__'))
        subprocess.run([sys.executable, '-B', str(target / 'verify.py')],
                       cwd=target, check=True)
        historical = json.loads((target / 'verification.json').read_text(encoding='utf-8'))
        assert (historical['points_36'], historical['points_37']) == (501633, 505359)
        target = Path(directory) / 'completion'
        shutil.copytree(ROOT / 'completion', target, ignore=shutil.ignore_patterns('__pycache__'))
        for script in ('verify.py', 'verify_coordinates.py', 'reproduce.py'):
            subprocess.run([sys.executable, '-B', str(target / script)], cwd=target, check=True)
        report = json.loads((target / 'verification.json').read_text(encoding='utf-8'))
        coordinates = json.loads((target / 'coordinate_verification.json').read_text(encoding='utf-8'))
        assert coordinates['all_checks_pass']
        for key in ('points_36', 'points_37', 'base_points', 'certificate_sha256'):
            assert coordinates[key] == report[key], key
    registered = results()
    assert report['all_checks_pass']
    assert report['unit_sphere_max_inner_product'] == '1/2'
    assert report['points_36'] == registered[0]['total']
    assert report['points_37'] == registered[1]['total']
    print(json.dumps(dict(status='PASS', results=registered,
                          construction='36D phase refill and 37D partial latitude orbits; two exact representations')))


if __name__ == '__main__':
    main()
