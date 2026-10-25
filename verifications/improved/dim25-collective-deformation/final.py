"""Derive the certificate count; run verify.py for the complete geometric proof."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    c=json.loads((ROOT/'source/collective_record25_certificate.json').read_text())
    circles=sum(len(k['owners']) for k in c['classes'])
    removed={tuple(u) for k in c['classes'] for u in k['owners']}
    removed.update(tuple(a+b for a,b in zip(u,k['shift'])) for k in c['classes'] for u in k['owners'])
    removed.update(map(tuple,c['local_owners']+c.get('new_owners',[])))
    total=196560-len(removed)+3*circles+len(c['small'])
    assert circles==971 and len(removed)==1991 and len(c['small'])==102
    assert total==c['points']==197584
    return [dict(dimension=25,total=total)]

if __name__=='__main__':
    for row in results(): print(json.dumps(row))
