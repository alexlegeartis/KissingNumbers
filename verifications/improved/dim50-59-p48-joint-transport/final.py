from pathlib import Path
import json
def results():
    rows = json.loads(Path(__file__).with_name('results.json').read_text())['rows']
    return [dict(dimension=r['dimension'],total=r['lower_bound']) for r in rows]
