"""Turn the discovered affine subcode into a small standalone certificate."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).parent;ROOT=HERE.parents[2]
SOURCE=ROOT/'verifications/improved/dim61-68-95-96-code-layers/source'
def expand(s,w):
    out=0;j=0
    for i in range(s.bit_length()):
        if s>>i&1:out|=((w>>j)&1)<<i;j+=1
    return out
def hamming(s,h,parity):
    positions=[i for i in range(s.bit_length()) if s>>i&1]
    order=[j for j,i in enumerate(positions) if h>>i&1]+[j for j,i in enumerate(positions) if not(h>>i&1)]
    assert len(positions)==16 and h.bit_count()==8 and h&s==h
    checks=[65535]+[sum(1<<order[i] for i in range(16) if i>>b&1) for b in range(4)]
    return dict(support=s,checks=checks,shift=(1<<order[0]) if parity else 0)
r=json.loads((HERE/'subcode67_best.json').read_text());context=json.loads((HERE/'subcode67_context.json').read_text())
phase=expand(r['support'],r['phase']);modified=[]
for kind,index,h in r['free']:
    s=context['old_order'][index] if kind==0 else context['extra_full'][index]
    modified.append(hamming(s,h,1^((phase&h).bit_count()%2)))
files=[SOURCE/'shared_parity/expanded67_68.npz',SOURCE/'dim61_96/code64_supports.npz',SOURCE/'dim61_96/cyclic73.weight16.bin',SOURCE/'dim61_96/shortening_support_search.json',SOURCE/'continuation49_72/netto.json']
out=dict(format='affine-sign-subcode-extension-v1',dimension=67,base_count=355090822,
    base_sign_family='shared_parity/expanded67_68.npz',source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    new_classes=[dict(support=r['support'],checks=r['checks'],shift=r['phase'],dimension=r['dimension'],count=r['count'])],
    modified_parent_classes=modified,total=355090822+r['count'])
(HERE/'subcode67_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(out['total'],len(modified))
