"""Derive totals from explicit class dimensions and inherited layer formulas.
Run verify.py (or --independent) for full affected-pair verification.
"""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
PREVIOUS={61:88975146,62:122539020,66:342271364}
def results():
    if not __debug__:raise RuntimeError('Run without -O.')
    totals={}
    for directory in ('dim61-62','dim66'):
        c=json.loads((HERE/directory/'certificate.json').read_text())
        f=c['family']; gain=sum(2**r['dimension'] for r in c['new_classes'])
        assert gain==c['gain'] and all(r['count']==2**r['dimension'] for r in c['new_classes'])
        n=len(f['new'])
        if directory=='dim61-62':
            assert f['removed_coordinates']==[24,42,53,55] and f['odd_total_parity'] and f['original_supports_retained']==9673
            totals[61]=(1<<26)+(9673+n)*2048+8555*16+122+gain
            totals[62]=6*(1<<24)+(9673+n)*2048+9145*16+124+gain
        else:
            assert not f['removed_coordinates'] and f['original_supports_retained']==30828
            totals[66]=(1<<28)+(34028+n)*2048+11088*16+132+gain
    return [dict(dimension=d,total=v,previous=PREVIOUS[d],gain=v-PREVIOUS[d]) for d,v in sorted(totals.items())]
if __name__=='__main__':print(json.dumps(results(),indent=2))
