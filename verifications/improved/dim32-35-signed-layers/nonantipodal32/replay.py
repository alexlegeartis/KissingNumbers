"""Portable isolated replay, complete manifest and non-antipodal controls."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: assertions are required.')
manifest = ROOT / 'SHA256SUMS.txt'
listed = set()
for line in manifest.read_text().splitlines():
    digest, name = line.split('  ', 1)
    path = (ROOT / name).resolve()
    assert ROOT in path.parents and path != manifest and name not in listed
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
    listed.add(name)
actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p != manifest}
assert listed == actual, 'manifest coverage'
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='nonantipodal32_replay_') as directory:
    target = Path(directory) / 'proof'
    shutil.copytree(ROOT, target, ignore=shutil.ignore_patterns('__pycache__'))
    for script in ('verify_coordinates.py', 'independent_audit.py', 'verify_fibers.py'):
        subprocess.run([sys.executable, '-B', str(target / script)], cwd=target, env=env, check=True)
    coordinate = json.loads((target / 'gain_coordinate_verification.json').read_text())
    report = json.loads((target / 'report.json').read_text())
    assert len(coordinate) == 2 and all(r['all_checks_pass'] and r['verified_points'] == 352992 for r in coordinate)
    assert report['status'] == 'PASS'
    assert {k: v['total'] for k, v in report['certificates'].items()} == {'fixed_body': 352992, 'discovery': 352992, 'baseline': 352832}
    assert report['delta']['net_gain'] == 160
    controls = report['negative_controls']
    assert len(controls) == 14 and all(r.get('rejected') or r.get('demonstrated') for r in controls)
    assert sum(r.get('rejected', False) for r in controls) == 13
    fixed = report['certificates']['fixed_body']
    assert fixed['materialized_distinct_rows'] == 352992 and fixed['materialized_global_unpaired_points'] == 380
    assert fixed['exceptional_pair_product_histogram']['-24'] == 768
    certificate = target / 'certificate32_352992_fixed_body_explicit.json'
    certificate.write_bytes(certificate.read_bytes() + b' ')
    bad = subprocess.run([sys.executable, '-B', str(target / 'independent_audit.py')], cwd=target, env=env, text=True, capture_output=True)
    assert bad.returncode != 0 and 'fixed_body SHA256' in bad.stderr, 'hash corruption must be rejected'
print(json.dumps(dict(status='PASS', dimension=32, total=352992, gain_over_published_repository=160, geometric_corruptions_rejected=13, absolute_product_false_rejection_demonstrated=True, hash_corruption_rejected=True)))
