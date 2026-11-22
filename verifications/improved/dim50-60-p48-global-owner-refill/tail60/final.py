import json
from pathlib import Path
def results():
    w=json.loads(Path(__file__).with_name('witness.json').read_text())
    return [dict(dimension=60,total=58636146+28*sum(q is not None for r in w['upgrade'] for q in r['tetrahedra']))]
