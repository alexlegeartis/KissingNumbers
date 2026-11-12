"""Derive current dimensions from explicit sign classes and exact section counts."""
import hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from geometry_counts import counts,verify_geometry
from section63 import verify_section

PREVIOUS={61:88260394,62:121824268,63:184286468,64:335496064,65:335564290,66:341994884,67:355090822}

def results():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    selection=json.loads((HERE/'data/selection.json').read_text())
    verify_geometry()
    dimensions=counts(HERE/selection['short'])
    dimensions.update(counts(HERE/selection['global']))
    dimensions[63]=verify_section(HERE/selection['global'],enumerate_words=False)['total']
    expanded=json.loads((HERE/selection['expanded66']).read_text())
    assert not expanded['removed_coordinates'] and expanded['original_supports_retained']==30828
    dimensions[66]=(1<<28)+(34028+len(expanded['new']))*2048+11088*16+132
    if 'subcode67' in selection:
        from subcode67 import results67
        dimensions[67]=results67(HERE/selection['subcode67'])
    return [dict(dimension=d,total=value,previous=PREVIOUS[d],gain=value-PREVIOUS[d]) for d,value in sorted(dimensions.items())]

if __name__=='__main__':print(json.dumps(results(),indent=2))
