#!/usr/bin/env python3
"""Look for an independently replayable clique partition proving the finite polar-pool ceiling."""
import json,sys,time
from pathlib import Path
import numpy as np
from baselines import DATA
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'compute17_20_20261104'))
from orbits17 import read_graph

off,nb=read_graph(DATA/'polar19_exact.graph');n=len(off)-1
adj=[set(map(int,nb[off[i]:off[i+1]])) for i in range(n)]
choices={v:sorted(u for u in adj[v] if u<1792) for v in range(1792,n)}
assert all(choices.values())
rng=np.random.default_rng(1912460);started=time.monotonic();success=None
for trial in range(400):
    groups=[[] for _ in range(1792)]
    order=sorted(range(1792,n),key=lambda v:(len(choices[v]),rng.random()))
    failed=False
    for v in order:
        valid=[u for u in choices[v] if all(w in adj[v] for w in groups[u])]
        if not valid:failed=True;break
        u=max(valid,key=lambda u:len(groups[u])+.2*rng.random())
        groups[u].append(v)
    if not failed:
        success=[[u]+g for u,g in enumerate(groups)];break
if success:
    assert sorted(v for g in success for v in g)==list(range(n))
    for g in success:
        assert all(w in adj[v] for j,v in enumerate(g) for w in g[j+1:])
    (DATA/'polar19_clique_partition.json').write_text(json.dumps(dict(groups=success),separators=(',',':'))+'\n')
report=dict(success=success is not None,trials=trial+1,seconds=time.monotonic()-started,
            vertices=n,cliques=len(success) if success else None)
(DATA/'polar19_cover_result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
