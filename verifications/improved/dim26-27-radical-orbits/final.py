"""Read the retained verified coordinate count."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
def results():
    with np.load(ROOT/'data/coordinates26_201182.npz',allow_pickle=False) as z:
        n,d=z['rational'].shape
        if (n,d)!=(201182,26):raise ValueError('Certificate dimensions')
    return [dict(dimension=d,total=n)]
if __name__=='__main__':
    for row in results():print(json.dumps(row))
