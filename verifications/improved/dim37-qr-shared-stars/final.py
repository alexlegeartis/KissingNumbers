"""Return the layered-code count for the registry."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent
def results():
    d=json.loads((P/'star_recipe.json').read_text())
    S=[(36,x) for x in d['partition']]+[(int(i),[*x,36]) for i,x in d['ordinary'].items()]
    E={tuple(sorted([i,j])) for i,x in S for j in x}
    n=533864+16*len(S)-4*len(E)
    assert n==534008
    n+=64*len(json.loads((P/'six_recipe.json').read_text())['ordinary_supports'])
    assert n==536440
    return [{'dimension':37,'total':n}]
if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
