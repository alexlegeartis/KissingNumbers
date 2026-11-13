from pathlib import Path
import numpy as np,json
ROOT=Path(__file__).resolve().parent

def results():
    rows=[]
    for d,n in [(29,231814),(30,248180)]:
        with np.load(ROOT/f'data/coordinates{d}_{n}.npz',allow_pickle=False) as z:
            if z['points'].shape!=(n,d):raise ValueError('Certificate dimensions')
        rows.append(dict(dimension=d,total=n))
    return rows
if __name__=='__main__':
    for r in results():print(json.dumps(r))
