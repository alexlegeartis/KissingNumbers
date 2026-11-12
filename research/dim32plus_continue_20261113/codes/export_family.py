"""Export explicit sign codes; search scores and pool indices are discarded."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from prepare_seed import SOURCE,DTYPE,shape
sys.path.insert(0,str(SOURCE/'shared_parity'))
from shared_layers import words,local,solve

HERE=Path(__file__).resolve().parent

def complete(s,required):
    basis=[65535]
    span={0,65535}
    for h in required:
        v=local(s,int(h))
        assert int(h)&s==int(h)
        if v in span:continue
        assert all((x^v).bit_count() in (0,8,16) for x in span)
        basis.append(v);span|={x^v for x in list(span)}
    assert len(basis)<=5
    rank=len(basis)
    columns=[sum(((v>>j)&1)<<k for k,v in enumerate(basis)) for j in range(16)]
    groups={v:[j for j in range(16) if columns[j]==v] for v in sorted(set(columns))}
    assert set(groups)==set(range(1,1<<rank,2))
    assert all(len(g)==1<<(5-rank) for g in groups.values())
    for bit in range(5-rank):
        h=sum(1<<j for group in groups.values() for i,j in enumerate(group) if i>>bit&1)
        basis.append(h)
    span=words(basis)
    assert len(set(span))==32 and {x.bit_count() for x in span}=={0,8,16}
    assert all(local(s,int(h)) in span for h in required)
    return basis

def export(prefix,poolname,removed,odd):
    old=list(map(int,np.fromfile(HERE/'old_supports.bin',dtype='<u8')))
    records=np.fromfile(HERE/poolname,dtype=DTYPE)
    raw=(HERE/f'{prefix}_best.txt').read_bytes()
    lines=raw.decode().splitlines();assert int(lines[0])==len(lines)-1
    added=[];constraints={}
    for line in lines[1:]:
        v=list(map(int,line.split()));assert len(v) in (7,8)
        s,phase,index=v[:3];pm=tuple(v[3:7]);extra=v[7] if len(v)==8 else 0
        assert int(records[index]['a'])==s and not s&removed
        H=shape(s,records[index]['old'],old)(pm)
        checks=complete(s,list(H)+([extra] if extra else []))
        shift=local(s,phase)
        assert not odd or shift.bit_count()%2==1
        added.append(dict(support=s,checks=checks,shift=shift))
        for j,h in zip(records[index]['old'],H):
            b=old[int(j)]
            if b&removed:continue
            assert h&(s&b)==h and h.bit_count()==8
            constraints.setdefault(b,[]).append((h,1^((h&phase).bit_count()%2)))
    changed=[]
    for s,equations in sorted(constraints.items()):
        checks=complete(s,[h for h,rhs in equations])
        eq=[(local(s,h),rhs) for h,rhs in equations]
        if odd:eq.append((65535,1))
        shift,_=solve(eq);assert shift is not None
        changed.append(dict(support=s,checks=checks,shift=shift))
    out=dict(format='explicit-hamming-sign-classes-v1',removed_coordinates=[i for i in range(64) if removed>>i&1],
             odd_total_parity=odd,new=sorted(added,key=lambda x:x['support']),old=changed,
             original_supports_retained=sum(not s&removed for s in old),
             source_archive_sha256=hashlib.sha256((SOURCE/'dim61_96/code64_supports.npz').read_bytes()).hexdigest(),
             discovery_seed_sha256=hashlib.sha256(raw).hexdigest())
    (HERE/f'{prefix}_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(certificate=f'{prefix}_certificate.json',new=len(added),modified_old=len(changed),
                         original=out['original_supports_retained'],odd=odd)),flush=True)

if __name__=='__main__':
    prefix=sys.argv[1]
    if prefix.startswith('short'):
        export(prefix,'short61_candidates.bin',sum(1<<c for c in (24,42,53,55)),True)
    elif prefix.startswith('expanded'):
        export(prefix,prefix+'_candidates.bin',0,False)
    else:
        export(prefix,'global64_candidates.bin',0,False)
