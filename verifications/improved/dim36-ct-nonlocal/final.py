"""Recompute the count for the exact 36D phase extension."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent
def results():
    c=json.loads((P/'certificate.json').read_text())
    n=len(c['outer_cosets'])-1429
    assert n==3671 and 403704+27*n==c['points']==502821
    return [{'dimension':36,'total':502821}]
if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
