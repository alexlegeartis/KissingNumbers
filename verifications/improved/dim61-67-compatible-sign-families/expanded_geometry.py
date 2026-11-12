"""Exact new-pair checks against the inherited 66D and 67D support layers."""
from pathlib import Path
import json,subprocess,tempfile
import numpy as np
from verify_family import SOURCE,rref,popcount

def check_expanded(certificate,dimension,cpp=True):
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    assert dimension in (66,67)
    cert=json.loads(Path(certificate).read_text());assert not cert['removed_coordinates']
    folder=SOURCE/'dim61_96'
    raw=np.fromfile(folder/'cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
    assert len(raw)==len(np.unique(raw,axis=0))==174324
    assert np.all(popcount(raw[:,0])+popcount(raw[:,1])==16)
    exps=(37,36,34,33,32,27,25,24,22,21,19,18,15,11,10,8,7,5,3,0)
    polynomial=sum(1<<i for i in exps);rows,pivots=rref([polynomial<<i for i in range(36)],73)
    for low,high in raw:
        word=int(low)|(int(high)<<64)
        for row,p in zip(rows,pivots):
            if word>>p&1:word^=row
        assert word==0
    choice=json.loads((folder/'shortening_support_search.json').read_text())['dimensions'][str(dimension)]
    retained=np.ones(len(raw),dtype=bool)
    for c in choice['deleted_coordinates']:retained&=((raw[:,c//64]>>np.uint64(c%64))&np.uint64(1))==0
    selected=raw[retained];assert len(selected)==choice['support_count']
    dual9=(21,37,38,49,57,63,67,70,72)
    head=[i for i in range(1,73) if i not in dual9]+[0]
    tails=[c for c in dual9 if c not in choice['deleted_coordinates']]
    assert len(head)==64 and len(tails)==dimension-64
    heads=np.zeros(len(selected),dtype=np.uint64)
    for j,c in enumerate(head):heads|=((selected[:,c//64]>>np.uint64(c%64))&np.uint64(1))<<np.uint64(j)
    ordinary=np.ones(len(selected),dtype=bool)
    for c in tails:ordinary&=((selected[:,c//64]>>np.uint64(c%64))&np.uint64(1))==0
    with np.load(folder/'code64_supports.npz',allow_pickle=False)as z:old=set(map(int,z['words']))
    assert len(heads[ordinary])==30828 and set(map(int,heads[ordinary]))==old
    extra=heads[~ordinary];assert len(extra)=={66:3200,67:9600}[dimension]
    coordinates=np.array([[(int(s)>>j)&1 for j in range(64)] for s in extra],dtype=np.int16)
    maximum=0
    for row in cert['new']:
        s=row['support'];packed=popcount(extra&np.uint64(s))
        physical=coordinates@np.array([(s>>j)&1 for j in range(64)],dtype=np.int16)
        assert np.array_equal(packed,physical) and int(packed.max())<=8
        maximum=max(maximum,int(packed.max()))
    netto=json.loads((SOURCE/'continuation49_72/netto.json').read_text())['dimensions'][str(dimension)]
    terms=dict(netto['terms']);assert terms['weight16']==len(selected)*2048
    terms['weight16']+=len(cert['new'])*2048
    previous={66:341994884,67:355090822}[dimension]
    scalar=None
    if cpp:
        with tempfile.TemporaryDirectory(prefix='expanded-sign-family-') as temporary:
            td=Path(temporary);executable=td/'check';data=td/'input.txt'
            lines=[f'{len(selected)} {len(cert["new"])}',' '.join(map(str,head))]
            lines += [f'{int(low)} {int(high)}' for low,high in selected]
            lines += [str(r['support']) for r in cert['new']]
            data.write_text('\n'.join(lines)+'\n')
            subprocess.run(['g++','-O3','-march=native','-std=c++17',str(Path(__file__).with_name('check_expanded.cpp')),'-o',str(executable)],check=True)
            result=subprocess.run([str(executable),str(data)],check=True,capture_output=True,text=True)
            scalar=json.loads(result.stdout)
            assert scalar['parent_support_pairs']==len(selected)*(len(selected)-1)//2
            assert scalar['new_extra_pairs']==len(extra)*len(cert['new'])
    return dict(status='PASS',dimension=dimension,total=sum(terms.values()),terms=terms,
                previous=previous,gain=sum(terms.values())-previous,
                parent_supports_checked=174324,retained_parent_supports=len(selected),extra_supports=len(extra),
                new_supports=len(cert['new']),extra_pairs=len(extra)*len(cert['new']),maximum_intersection=maximum,
                scalar_cpp=scalar,
                input_scope='Every retained parent-support pair is checked directly, so the inherited 73-bit code minimum is unnecessary. The unchanged Netto four-support construction is replayed by verify.py.')

if __name__=='__main__':
    import sys
    print(json.dumps(check_expanded(sys.argv[1],int(sys.argv[2])),indent=2))
