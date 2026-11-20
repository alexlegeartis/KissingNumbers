from pathlib import Path
import json

def results():
    p=Path(__file__).resolve().parent
    d=json.loads((p/'verification.json').read_text())
    c=json.loads((p/'certificate.json').read_text())
    assert d['status']=='PASS' and d['points']==c['claimed_size']==774876
    assert d['dense']+d['octad_points']+d['roots']+d['stars']==d['points']
    return [dict(dimension=39,total=d['points'])]
