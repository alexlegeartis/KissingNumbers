"""Prepare the fixed published 67D family without changing its files."""
from pathlib import Path
import json,numpy as np
from review_codes import expand,words,HERE,CODES,SOURCE
with np.load(SOURCE/'shared_parity/expanded67_68.npz',allow_pickle=False) as z:
    data={}
    for kind in ('old','new'):
        for s,h,t in zip(z[kind+'_supports'],z[kind+'_duals'],z[kind+'_shifts']):
            s=int(s);h=list(map(int,h));data[s]=(expand(s,int(t)),[expand(s,w) for w in words(h) if w.bit_count()==8])
    selected=list(map(int,z['new_supports']))
old=list(map(int,np.fromfile(CODES/'old_supports.bin',dtype='<u8')))
with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False) as z:assert set(old)==set(map(int,z['words']))
raw=np.fromfile(SOURCE/'dim61_96/cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
choice=json.loads((SOURCE/'dim61_96/shortening_support_search.json').read_text())['dimensions']['67'];deleted=choice['deleted_coordinates']
dual9=(21,37,38,49,57,63,67,70,72);head=[i for i in range(1,73) if i not in dual9]+[0];tail=[i for i in dual9 if i not in deleted]
coords=head+tail;assert len(coords)==67
allwords=[int(lo)|(int(hi)<<64) for lo,hi in raw];kept=[w for w in allwords if all(not(w>>i&1) for i in deleted)]
full=[sum(((w>>c)&1)<<j for j,c in enumerate(coords)) for w in kept];extra=[w for w in full if w>>64]
assert len(extra)==9600 and len(kept)==40428
lines=[f'{len(old)} {len(selected)} {len(extra)}']
for s in old+selected:
    phase,checks=data.get(s,(0,[]));lines.append(' '.join(map(str,[s,phase,len(checks),*checks])))
lines.append(' '.join(str(w&((1<<64)-1)) for w in extra))
(HERE/'subcode67_input.txt').write_text('\n'.join(lines)+'\n')
(HERE/'subcode67_context.json').write_text(json.dumps(dict(old_order=old,extra_full=extra,coordinates=coords,published_new=selected)))
print('prepared',len(old),len(selected),len(extra))
