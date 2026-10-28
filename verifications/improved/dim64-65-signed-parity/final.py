"""Read the explicit support family and derive both construction sizes."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def results():
    data = json.loads((ROOT/'data/construction.json').read_text())
    base = 2**28 + 30828*2048 + 10416*16 + 128
    n64 = base + 2048*len(data['added_blocks'])
    n65 = n64 + 65536 + 336*8 + 2
    assert data['dimension64'] == n64 == 332170112
    assert data['dimension65'] == n65 == 332238338
    path = ROOT.parent/'dim61-68-95-96-code-layers/source/shared_parity/shared_layers.json'
    shared = json.loads(path.read_text())
    assert shared['status'].startswith('Verified:')
    current = shared['dimensions']
    n64 = base+2048*shared['families']['global']['new_supports']
    n65 = n64+data['owner_points65']+2
    assert current['64']['total'] == n64 and current['65']['total'] == n65
    assert current['64']['gain'] == n64-data['dimension64'] > 0
    assert current['65']['gain'] == n65-data['dimension65'] > 0
    return [dict(dimension=64,total=n64),dict(dimension=65,total=n65)]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
