"""Derive the exact certificate count; run verify.py for geometry."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

def results():
    cert=json.loads((ROOT/'source/certificate.json').read_text())
    stars=cert['heavy_stars']
    pairs={tuple(sorted((s['center'],b))) for s in stars for b in s['leaves']}
    kernel=[s for s in (ROOT/'source/source/dense_kernel.txt').read_text().splitlines() if s and s[0] not in '#$']
    reps=[s for s in (ROOT/'source/source/dense_representatives.txt').read_text().splitlines() if s and s[0] not in '#$']
    supports=(ROOT/'source/source/supports_3424.txt').read_text().splitlines()
    total=len(reps)*(1<<len(kernel))+128*len(supports)+4*(39*38//2-len(pairs))+16*len(stars)
    assert total==cert['claimed_size'] and cert['dimension']==39
    return [dict(dimension=39,total=total)]

if __name__ == '__main__':
    for row in results(): print(json.dumps(row))
