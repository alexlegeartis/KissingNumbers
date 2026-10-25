"""Integrity checks and isolated exact replay of both midpoint certificates."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, temporary_source, replay
from final import results


if __name__ == '__main__':
    check_manifest(ROOT)
    expected = results()
    with temporary_source(ROOT) as source:
        for row in expected:
            record = json.loads((source / ('record%d.json' % row['dimension'])).read_text())
            replay(source, 'verify.py', record['filename'])
    print(json.dumps(dict(status='PASS', results=expected)))
