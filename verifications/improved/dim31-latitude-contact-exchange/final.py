"""Count the supplied configuration; verify.py establishes its geometry."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    row=json.loads((ROOT/'data/provenance.json').read_text())['records'][0]
    x=np.load(ROOT/row['certificate'],allow_pickle=False)['points']
    if x.shape!=(row['total'],31) or row['total']!=row['previous_repository']+row['gain']:
        raise ValueError('Count/dimension/gain mismatch')
    return [{'dimension':31,'total':len(x)}]
if __name__=='__main__':
    for row in results():print(json.dumps(row))
