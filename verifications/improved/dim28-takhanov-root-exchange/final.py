"""Read the frozen coordinate certificate used by the repository registry."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
def results():
    m=json.loads((ROOT/'manifest.json').read_text())
    with np.load(ROOT/m['coordinate_file'],allow_pickle=False)as z:
        n,d=z['points'].shape
    if (n,d)!=(m['count'],m['dimension']):raise ValueError('frozen coordinate dimensions')
    return [dict(dimension=d,total=n)]
if __name__=='__main__':
    for row in results():print(json.dumps(row))
