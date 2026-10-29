"""Complete geometry, direct-integer insertion checks and integrity replay."""
from pathlib import Path
import json
import sys
from verify_geometry import verify
from verify_delta import replay

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest


if __name__ == '__main__':
    check_manifest(ROOT)
    wanted = set(map(int, sys.argv[1:]))
    for record in json.loads((ROOT/'data/provenance.json').read_text())['records']:
        if wanted and record['dimension'] not in wanted:
            continue
        actual = verify(ROOT/record['certificate'])
        expected = json.loads((ROOT/record['verification']).read_text())
        assert {k:v for k,v in actual.items() if k != 'seconds'} == {
            k:v for k,v in expected.items() if k != 'seconds'}
        delta = replay(record)
        assert delta == record['direct_integer_replay']
        print(json.dumps(dict(status='PASS', geometry=actual, delta=delta)), flush=True)
