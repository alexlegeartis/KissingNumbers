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
    return [dict(dimension=64,total=n64),dict(dimension=65,total=n65)]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
