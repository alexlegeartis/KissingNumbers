"""Verify complete integrity, then replay all proofs in a temporary directory."""
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
    raise SystemExit('Run without -O: assertions must be enabled.')
manifest = ROOT / 'SHA256SUMS.txt'
listed = set()
for line in manifest.read_text().splitlines():
    digest, name = line.split('  ', 1)
    path = (ROOT / name).resolve()
    assert ROOT in path.parents and path != manifest and name not in listed
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
    listed.add(name)
actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
          if p.is_file() and '__pycache__' not in p.parts and p != manifest}
assert listed == actual, 'manifest coverage'
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='nonantipodal33_replay_') as directory:
    target = Path(directory) / 'proof'
    shutil.copytree(ROOT, target, ignore=shutil.ignore_patterns('__pycache__'))
    producer = []
    for script in ('verify.py', 'verify_coordinates.py'):
        p = subprocess.run([sys.executable, '-B', str(target/script)], cwd=target,
                           env=env, text=True, capture_output=True, check=True)
        result = json.loads(p.stdout)
        assert result['all_checks_pass'] and result['points'] == 370408
        assert result['certificate_sha256'] == '5686dfcd2b891bf6d9be6b01013d8205774aa1a19c02a2c45bdfd86194314a60'
        producer.append(result)
    assert producer[0]['unpaired_points'] == 500
    assert producer[1]['rank_mod_101'] == 33
    for script, output in [('independent_core.py', 'intermediate'),
                           ('independent_audit.py', 'final')]:
        subprocess.run([sys.executable, '-B', str(target/script), '--output',
                        str(Path(directory)/output)], cwd=target, env=env, check=True)
    middle = json.loads((Path(directory)/'intermediate/report.json').read_text())
    final = json.loads((Path(directory)/'final/report.json').read_text())
    for report, count in [(middle, 370392), (final, 370408)]:
        assert report['status'] == 'PASS' and report['lower_bound'] == count
        assert report['candidate']['total'] == count
        assert report['candidate']['rank_mod101'] == 33
        assert report['controls']['rejection_count'] == 16
        assert len(report['controls']['results']) == 16
        assert all(r['rejected'] for r in report['controls']['results'])
    assert middle['old_control']['total'] == 370376
    assert middle['delta']['net_gain'] == 16
    assert final['delta']['gain_vs_370392'] == 16
    assert final['delta']['gain_vs_actual_370376'] == 32
    assert final['delta']['actual_old_control_added'] == 1796
    assert final['delta']['actual_old_control_removed'] == 1764
    assert final['candidate']['counts'] == {'dense':131072, 'full':235136,
           'explicit':2064, 'root':1956, 'star':180}
    assert final['candidate']['missing_antipodes'] == 500
    assert final['delta']['missing_antipodes_by_category'] == {'explicit':372, 'root':12, 'star':116}
    # Retained reports are evidence only; independently regenerated row digests
    # must also agree with the frozen audit's complete coordinate materialization.
    for report, name in [(middle, 'original_intermediate_report.json'),
                         (final, 'original_final_report.json')]:
        archived = json.loads((target/'evidence'/name).read_text())
        assert report['candidate']['row_sha256'] == archived['candidate']['row_sha256']
    certificate = target/'certificate33_370408_explicit.json'
    certificate.write_bytes(certificate.read_bytes() + b' ')
    p = subprocess.run([sys.executable, '-B', str(target/'independent_audit.py'),
                        '--output', str(Path(directory)/'corrupt')], cwd=target,
                       env=env, text=True, capture_output=True)
    assert p.returncode != 0 and 'SHA-256 mismatch' in p.stderr, 'hash mutation accepted'
print(json.dumps(dict(status='PASS', dimension=33, total=370408,
      gain_over_actual_published_repository=32, intermediate_not_previously_published=370392,
      missing_antipodes=500, intermediate_controls=16, final_controls=16,
      additional_hash_corruption_rejected=True)))
