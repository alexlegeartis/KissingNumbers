"""Derive the current certificate count; verify.py proves the geometry."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

def results():
    c=json.loads((ROOT/'refill/certificate.json').read_text())
    assert c['format']=='fano_grid_cube_latitudes_v1'
    assert len(c['body_line_indices'])==len(c['doubled_head_lines'])==1568
    assert len(c['removed_cap_body_lines'])==20
    assert len(c['additional_latitude_points'])==384
    total=589680-6*len(c['removed_cap_body_lines'])+13710+len(c['axis_keep'])+2*len(c['additional_coordinate_axes'])+4*len(c['doubled_head_lines'])+len(c['additional_latitude_points'])
    assert c['dimension']==38 and total==c['claimed_size']==611522
    return [dict(dimension=38,total=total)]

if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
