"""Derive the certificate count; run verify.py for the complete geometric proof."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent

def results():
    c=json.loads((ROOT/'source/certificate.json').read_text())
    total=589680+13710+len(c['axis_keep'])+2*len(c['additional_coordinate_axes'])+4*len(c['doubled_head_lines'])
    assert len(c['body_line_indices'])==len(c['doubled_head_lines'])==1352
    assert c['dimension']==38 and total==c['claimed_size']==610394
    return [dict(dimension=38,total=total)]

if __name__=='__main__':
    for row in results(): print(json.dumps(row))
