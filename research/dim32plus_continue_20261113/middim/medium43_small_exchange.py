"""Exact one-for-two and two-for-three searches on a prepared finite pool."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import json,sys,time
import numpy as np
HERE=Path(__file__).parent
def run(extra=False,maxremoved=2):
    stem='medium43_pool_f59'+('_extra' if extra else '')
    P=json.loads((HERE/f'{stem}.json').read_text());F=np.load(HERE/f'{stem}.npz')
    B=F['selfbeta'];K=F['paircut'];V=P['vertices'];seen=set();groups=defaultdict(list)
    for i,v in enumerate(V):
        key=(v['row'],tuple(v['blockers']),tuple(K[v['head']]))
        if key in seen:continue
        seen.add(key)
        if len(v['blockers'])<=maxremoved:groups[tuple(v['blockers'])].append(i)
    def pair(i,j):return B[V[i]['row'],V[j]['row']]<=K[V[i]['head'],V[j]['head']]
    def emit(ids,deleted):
        result=dict(status='Exact finite exchange; independent geometric replay required',source_pool=stem,gain=len(ids)-len(deleted),new_vertex_indices=ids,deleted_current_indices=deleted)
        (HERE/(stem+'_small_exchange.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True);return result
    if () in groups:return emit(groups[()][:1],[])
    for deleted,ids in groups.items():
        if len(deleted)!=1:continue
        for a,b in combinations(ids,2):
            if pair(a,b):return emit([a,b],list(deleted))
    tries=0;start=time.time()
    for deleted,ids in groups.items():
        if len(deleted)!=2:continue
        choices=ids+groups.get((deleted[0],),[])+groups.get((deleted[1],),[])
        choices=list(dict.fromkeys(choices))
        adjacency=[set() for _ in choices]
        for i,a in enumerate(choices):
            for j in range(i):
                if pair(a,choices[j]):adjacency[i].add(j);adjacency[j].add(i)
        for i in range(len(choices)):
            for j in adjacency[i]:
                if j<=i:continue
                for k in adjacency[i]&adjacency[j]:
                    if k>j:return emit([choices[i],choices[j],choices[k]],list(deleted))
        tries+=1
    print(json.dumps(dict(status='No gain in this finite one/two-deletion test',groups=len(groups),two_groups=tries,seconds=time.time()-start)),flush=True)
if __name__=='__main__':run('--extra' in sys.argv)
