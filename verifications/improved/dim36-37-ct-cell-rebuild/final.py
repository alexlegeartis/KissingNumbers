"""Return the count of the retained complete-cell construction."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent

def results():
    c=json.loads((P/'certificate.json').read_text())
    assert c['dimension']==36 and c['squared_norm']==12
    n=3*756+756*(3*81)+672*(18*18)+(len(c['outer_cosets'])-1429)*27
    assert n==c['points']==502767
    n37=n+18*len(c['lift_cosets'])
    assert n37==c['points_37']
    return [{'dimension':36,'total':n},{'dimension':37,'total':n37}]

if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
