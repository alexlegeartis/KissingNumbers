"""Regenerate exact two-copy split options from the prepared proposal pool."""
from pathlib import Path
import json,sys,numpy as np
HERE=Path(__file__).parent
extra='--extra' in sys.argv;suffix='_extra' if extra else '';limit=4 if extra else 100
p=json.loads((HERE/f'medium43_pool_f59{suffix}.json').read_text());f=np.load(HERE/f'medium43_pool_f59{suffix}.npz')
V=p['vertices'];K=f['paircut'];seen={};byrow={}
for i,v in enumerate(V):
    if len(v['blockers'])>limit:continue
    key=(v['row'],tuple(v['blockers']),tuple(K[v['head']]))
    if key in seen:continue
    seen[key]=i;byrow.setdefault(v['row'],[]).append(i)
options=[];signatures=set()
for row,ids in byrow.items():
    for a,i in enumerate(ids):
        for j in ids[:a]:
            if K[V[i]['head'],V[j]['head']]<24:continue
            blocked=tuple(sorted(set(V[i]['blockers'])|set(V[j]['blockers'])))
            if len(blocked)>limit:continue
            signature=(row,blocked,tuple(np.minimum(K[V[i]['head']],K[V[j]['head']])))
            if extra and signature in signatures:continue
            signatures.add(signature);options.append([i,j,list(blocked)])
(HERE/f'medium43_options{suffix}.json').write_text(json.dumps(dict(options=options,deduplicated_vertices=list(seen.values()))))
print('deduplicated proposals',len(seen),'options',len(options))
