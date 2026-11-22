import json
from pathlib import Path

def results():
    """Registry-facing values; run verify.py for constructive replay."""
    x=json.loads(Path(__file__).with_name('results.json').read_text())
    return [dict(dimension=int(d),total=t) for d,t in sorted(x.items(),key=lambda i:int(i[0]))]
