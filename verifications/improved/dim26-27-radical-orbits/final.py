"""Read the retained verified coordinate count."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
def results():
    rows=[]
    for dimension,count in ((26,201188),(27,209884)):
        with np.load(ROOT/f'data/coordinates{dimension}_{count}.npz',allow_pickle=False) as z:
            n,d=z['rational'].shape
            if (n,d)!=(count,dimension):raise ValueError('Certificate dimensions')
        rows.append(dict(dimension=d,total=n))
    return rows
if __name__=='__main__':
    for row in results():print(json.dumps(row))
