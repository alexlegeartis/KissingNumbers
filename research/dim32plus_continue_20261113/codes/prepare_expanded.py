"""Filter explicit candidates against the inherited 66/67D extra supports."""
import itertools,json,subprocess
from pathlib import Path
import numpy as np
from prepare_seed import SOURCE,HERE,DTYPE,shape
from verify_family import span,expand

def extra_heads(d):
    mother=np.fromfile(SOURCE/'dim61_96/cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
    deleted=(21,37,38,49,57,63,67,70,72)
    head=[i for i in range(1,73) if i not in deleted]+[0]
    choice=json.loads((SOURCE/'dim61_96/shortening_support_search.json').read_text())['dimensions'][str(d)]
    kept=np.ones(len(mother),dtype=bool)
    for c in choice['deleted_coordinates']:kept&=((mother[:,c//64]>>np.uint64(c%64))&np.uint64(1))==0
    mother=mother[kept];assert len(mother)==choice['support_count']
    heads=np.zeros(len(mother),dtype=np.uint64)
    for j,c in enumerate(head):heads|=((mother[:,c//64]>>np.uint64(c%64))&np.uint64(1))<<np.uint64(j)
    with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False)as z:old=set(map(int,z['words']))
    assert sum(int(s) in old for s in heads)==30828
    return np.array([s for s in heads if int(s) not in old],dtype='<u8')

def seed(d):
    extra=extra_heads(d);extra.tofile(HERE/f'extra{d}.bin')
    subprocess.run([str(HERE/'filter_expanded.bin'),str(HERE/'candidates.bin'),str(HERE/f'extra{d}.bin'),str(HERE/f'expanded{d}_candidates.bin')],check=True)
    pool=np.fromfile(HERE/f'expanded{d}_candidates.bin',dtype=DTYPE);ids={int(r['a']):i for i,r in enumerate(pool)}
    old=np.fromfile(HERE/'old_supports.bin',dtype='<u8')
    global_rows=json.loads((HERE/'global_repaired_certificate.json').read_text())['new']
    global_rows=[r for r in global_rows if r['support'] in ids]
    name='expanded66' if d==66 else 'expanded67_68'
    with np.load(SOURCE/f'shared_parity/{name}.npz',allow_pickle=False)as z:
        published=[dict(support=int(s),checks=list(map(int,h)),shift=int(t)) for s,h,t in zip(z['new_supports'],z['new_duals'],z['new_shifts'])]
    chosen=global_rows if len(global_rows)>=len(published) else published
    lines=[]
    for r in chosen:
        s=r['support'];i=ids[s];dual={expand(s,h) for h in span(r['checks'])};check=shape(s,pool[i]['old'],old)
        pm=next(pm for pm in itertools.permutations(range(4)) if set(check(pm))<=dual)
        partial={0,s}
        for h in check(pm):partial|={x^h for x in list(partial)}
        assert len(partial)==16
        completion=min(dual-partial)
        lines.append(' '.join(map(str,(s,expand(s,r['shift']),i,*pm,completion))))
    (HERE/f'expanded{d}_best.txt').write_text(str(len(lines))+'\n'+'\n'.join(lines)+'\n')
    return dict(dimension=d,candidates=len(pool),extra_supports=len(extra),global_seed=len(global_rows),published_seed=len(published),chosen=len(chosen))

if __name__=='__main__':
    subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'filter_expanded.cpp'),'-o',str(HERE/'filter_expanded.bin')],check=True)
    result=[seed(d) for d in (66,67)]
    (HERE/'expanded_preparation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
