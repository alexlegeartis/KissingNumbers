"""Recompute layer cardinalities and explicit 63D section, without receipts.
Full geometric verification is provided by verify.py and audit/.
"""
from pathlib import Path
from itertools import combinations
from math import comb
import json
import numpy as np
HERE=Path(__file__).resolve().parent
PREVIOUS={63:184393988,64:335711104,65:335779330}

def _words(checks):
    rows={}
    for h in checks:
        x=h
        while x:
            j=x.bit_length()-1
            if j in rows:x ^= rows[j]
            else:rows[j]=x;break
    basis=[]
    for j in range(16):
        if j in rows:continue
        w=1<<j
        for k in sorted(rows):
            if (w&rows[k]).bit_count()%2:w ^= 1<<k
        assert all((w&h).bit_count()%2==0 for h in checks)
        basis.append(w)
    out=[0]
    for w in basis:out += [v^w for v in out]
    return np.asarray(out,dtype=np.uint16)

def results():
    if not __debug__:raise RuntimeError('Run without -O.')
    p=HERE.parent/'dim61-67-compatible-sign-families'
    source=HERE.parent/'dim61-68-95-96-code-layers/source'
    c=json.loads((HERE/'certificate.json').read_text())
    f=json.loads((p/'data/global64.json').read_text())
    with np.load(source/'dim61_96/code64_supports.npz',allow_pickle=False) as z:
        old=list(map(int,z['words']))
        gen=[sum(int(v)<<i for i,v in enumerate(r)) for r in z['generator']]
    assert len(old)==30828 and len(f['new'])==1940
    rows={r['support']:r for r in f['old']+f['new']+c['modified_parent_classes']+c['new_classes']}
    words={s:_words(r['checks'])^r['shift'] for s,r in rows.items()}
    gain=sum(len(words[r['support']]) for r in c['new_classes'])
    base=(1<<28)+(len(old)+len(f['new']))*2048+10416*16+128
    total64=base+gain
    # Construction choice, not a cached cardinality.
    h=9606450284331136
    assert h in json.loads((source/'shared_parity/dual_octads.json').read_text())
    pos=[i for i in range(64) if h>>i&1]
    assert len(pos)==8 and all((g&h).bit_count()%2==0 for g in gen)
    projection=[sum(((g>>i)&1)<<j for j,i in enumerate(pos)) for g in gen]
    piv={}
    for x in projection:
        while x:
            j=x.bit_length()-1
            if j in piv:x ^= piv[j]
            else:piv[j]=x;break
    assert len(piv)==7
    def count(s):
        m=(s&h).bit_count()
        if m%2:return 0
        if s not in rows:
            assert m<=8
            return 1120 if m==8 else comb(m,m//2)*(1<<(11-m))
        mask=sum(((h>>i)&1)<<j for j,i in enumerate(i for i in range(64) if s>>i&1))
        return int(np.count_nonzero(np.bitwise_count(words[s]&mask)==m//2))
    blocks={tuple(sorted((a,b,b2,a^b^b2))) for a,b,b2 in combinations(range(64),3)}
    assert len(blocks)==10416
    four=sum(comb(m,m//2)*(1<<(4-m)) for t in blocks if (m:=sum(i<8 for i in t))%2==0)
    total63=comb(8,4)*(1<<(28-7))+sum(count(s) for s in old)+sum(count(r['support']) for r in f['new']+c['new_classes'])+four+2*(64-8)
    totals={63:total63,64:total64,65:total64+68224+2}
    return [dict(dimension=d,total=v,previous=PREVIOUS[d],gain=v-PREVIOUS[d]) for d,v in sorted(totals.items())]

if __name__=='__main__':print(json.dumps(results(),indent=2))
