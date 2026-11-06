"""Derive the current count; verify.py proves the complete geometry."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
def results():
 c=json.loads((ROOT/'hamming/certificate.json').read_text());old=json.loads((ROOT/'refill/certificate.json').read_text());assert c['format']=='orthogonal_hamming_38_v1' and len(c['removed_cap_body_lines'])==22 and c['new_point_count']==512
 total=6*(98280-len(c['removed_cap_body_lines']))+13710+len(old['axis_keep'])+2*len(old['additional_coordinate_axes'])+4*len(old['body_line_indices'])+c['new_point_count'];assert c['dimension']==38 and total==c['claimed_size']==611638;return [dict(dimension=38,total=total)]
if __name__=='__main__':
 for row in results():print(json.dumps(row))
