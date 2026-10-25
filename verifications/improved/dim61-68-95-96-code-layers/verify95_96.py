"""Fast independent replay of the current 95--96D support certificate."""
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results

check_manifest(ROOT)
with temporary_source(ROOT) as source:
    replay(source, 'dim61_96/mixed_pair_certificate.py')
    fresh = json.loads((source/'dim61_96/mixed_pair_certificate.json').read_text())
    expected = {r['dimension']:r['total'] for r in results() if r['dimension'] in (95, 96)}
    assert {int(d):r['bound'] for d,r in fresh['dimensions'].items()} == expected
print(json.dumps({'status':'PASS','results':expected}))
