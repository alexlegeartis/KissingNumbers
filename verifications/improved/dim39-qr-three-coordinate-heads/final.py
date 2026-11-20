from pathlib import Path
import json
def results():
    p=Path(__file__).resolve().parent
    r=json.loads((p/'verification.json').read_text())
    s=json.loads((p/'certificate.json').read_text())
    assert r['status']=='PASS' and r['points']==s['claimed_size']==783468
    assert r['points']==sum(r[k] for k in ('dense','octad_points','roots','stars'))
    return [dict(dimension=39,total=r['points'])]
