"""Derive the certificate count; run verify.py for the complete geometric proof."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    c=json.loads((ROOT/'source/orbit_certificate.json').read_text())
    b=json.loads((ROOT/'source/base_certificate.json').read_text())
    stars=b['heavy_stars']
    pairs={tuple(sorted((s['center'],j))) for s in stars for j in s['leaves']}
    supports=(ROOT/'source/source/supports_3424.txt').read_text().splitlines()
    sparse=128*len(supports)+4*(39*38//2-len(pairs))+16*len(stars)
    assert sparse==441244 and len(c['group_words'])==256 and len(c['representatives'])==1281
    total=sparse+len(c['group_words'])*len(c['representatives'])
    assert c['dimension']==39 and total==c['claimed_size']==769180
    return [dict(dimension=39,total=total)]

if __name__=='__main__':
    for row in results(): print(json.dumps(row))
