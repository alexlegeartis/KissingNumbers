"""Conservative published lower bounds; see README.md and verify.py."""
import json
from pathlib import Path
def results():
    return [dict(dimension=int(d),total=n) for d,n in json.loads((Path(__file__).parent/'claims.json').read_text()).items()]
if __name__=='__main__':
    for row in results():print(json.dumps(row))
