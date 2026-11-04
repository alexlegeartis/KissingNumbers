"""Count the supplied explicit configurations; verify.py proves their geometry."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    rows=json.loads((ROOT/'data/provenance.json').read_text())['records']
    result=[]
    for row in rows:
        points=np.load(ROOT/row['certificate'],allow_pickle=False)['points']
        if points.shape!=(row['total'],row['dimension']):raise ValueError('Coordinate count/dimension')
        if row['total']!=row['previous_repository']+row['gain']:raise ValueError('Gain arithmetic')
        result.append({'dimension':row['dimension'],'total':len(points)})
    return result

if __name__=='__main__':
    for row in results():print(json.dumps(row))
