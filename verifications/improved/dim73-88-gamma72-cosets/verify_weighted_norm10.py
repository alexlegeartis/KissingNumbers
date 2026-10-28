"""Fast isolated replay of the mixed-tail weighted-priority existence proofs."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results

check_manifest(ROOT)
with temporary_source(ROOT) as source:
    replay(source, 'norm10/weighted_norm10_certificate.py', '--check')
    replay(source, 'norm10/verify_weighted_independent.py')
    report = json.loads((source/'norm10/weighted_norm10_certificate.json').read_text())
    fresh = {int(d): r['bound'] for d, r in report['dimensions'].items()}
    expected = {r['dimension']: r['total'] for r in results() if r['dimension'] in fresh}
    assert fresh == expected
print(json.dumps({'status': 'PASS', 'results': fresh}))
