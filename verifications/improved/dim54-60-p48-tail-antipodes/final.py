from pathlib import Path
import json
def results():
    rows = json.loads(Path(__file__).with_name('claims.json').read_text())
    return [dict(dimension=int(d),total=v['total']) for d,v in sorted(rows.items())]
