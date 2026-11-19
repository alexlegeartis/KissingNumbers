"""Count the triangle layer and inherited layers for the results registry."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent

def results():
    parent=P.parent/'dim38-joint-latitudes'
    old=json.loads((parent/'refill/certificate.json').read_text())
    h=json.loads((parent/'hamming/certificate.json').read_text())
    assert h['claimed_size']==611638 and h['new_point_count']==512
    n=6*(98280-14-24)+2*6855+len(old['axis_keep'])+2*len(old['additional_coordinate_axes'])+4*len(old['body_line_indices'])+2*128*3
    assert n==611798
    return [{'dimension':38,'total':n}]

if __name__=='__main__':
    for row in results():
        print(json.dumps(row))
