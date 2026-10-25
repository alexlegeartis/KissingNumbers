"""Derive the exact certificate count; run verify.py for geometry."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

def results():
    cert=json.loads((ROOT/'source/certificate.json').read_text())
    new=4*len(cert['doubled_head_lines'])
    assert len(cert['body_line_indices'])==new//4
    total=589680+13710+len(cert['axis_keep'])+new
    assert total==cert['claimed_size'] and cert['dimension']==38
    return [dict(dimension=38,total=total)]

if __name__ == '__main__':
    for row in results(): print(json.dumps(row))
