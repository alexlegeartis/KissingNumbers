"""Report the size of the verified construction."""
import json
from pathlib import Path
def results():
    v = json.loads(Path(__file__).with_name('verification.json').read_text())
    assert v['status'] == 'PASS' and v['dimension'] == 38
    return [{'dimension':38, 'total':v['points']}]
if __name__ == '__main__':
    for x in results():
        print(json.dumps(x))
