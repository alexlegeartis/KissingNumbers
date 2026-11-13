"""Bounded exact E8/D7 orientations and optimal two-copy unions for LS15."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
spec=importlib.util.spec_from_file_location('_bw_orientation',REPO/'research/reclaim_26_31_20261112/modules/bw15_transplant.py')
oldmodule=importlib.util.module_from_spec(spec);spec.loader.exec_module(oldmodule)
NORMAL=oldmodule.NORMAL
EMBED=np.r_[NORMAL,np.arange(24,31)]


def extract(graph,match):
    reverse=np.full(graph.shape[1],-1,np.int32)
    for i,j in enumerate(match):
        if j>=0:reverse[j]=i
    left=set(map(int,np.flatnonzero(match<0)));right=set();queue=list(left)
    for i in queue:
        for j0 in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
            j=int(j0)
            if j==match[i] or j in right:continue
            right.add(j)
            if reverse[j]>=0 and int(reverse[j]) not in left:
                left.add(int(reverse[j]));queue.append(int(reverse[j]))
    return sorted(left),sorted(right)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path);ap.add_argument('module',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=90);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();x=np.load(a.source,allow_pickle=False)['points'].astype(np.int16)//2
    m=np.load(a.module,allow_pickle=False)['points'].astype(np.int16);base=m[np.sum(m[:,8:].astype(np.int32)**2,axis=1)==16]
    assert base.shape==(2240,15)
    q=np.sum(x[:,24:].astype(np.int32)**2,axis=1);h=np.sum(x[:,NORMAL].astype(np.int32)**2,axis=1)
    relevant=((q==8)&(h>8))|((q==16)&(h>0));ids=np.flatnonzero(relevant)
    old=x[ids][:,EMBED];oldf=old.astype(np.float32)
    outside=[i for i in range(24) if i not in NORMAL]
    external=np.any(x[ids][:,outside],axis=1)
    allmodule=~np.any(x[:,outside],axis=1)
    oldset={r.tobytes() for r in x[allmodule][:,EMBED]}
    rng=np.random.default_rng(20261113);best=-1;records=[];trial=0
    while trial==0 or time.monotonic()-start<a.seconds:
        y,orientation=oldmodule.orientation(base,trial,rng)
        new=np.array([r for r in y if r.tobytes() not in oldset],dtype=np.int16)
        g=csr_matrix(new.astype(np.float32)@oldf.T>16)
        match=maximum_bipartite_matching(g,perm_type='column')
        gain=len(new)-int(np.sum(match>=0))
        cost=int(np.sum(np.any(y.astype(np.float32)@oldf[external].T>16,axis=0)))
        record=dict(trial=trial,gain=gain,new_points=len(new),edges=g.nnz,external_cost=cost,whole_module_gain=282-cost)
        if gain>best:
            added,removed=extract(g,match);remove=ids[removed];add=np.zeros((len(added),31),np.int16);add[:,EMBED]=new[added]
            keep=np.ones(len(x),bool);keep[remove]=False;f=np.r_[x[keep],add]
            assert len(f)==len(x)+gain and len(f)==len(np.unique(f,axis=0))
            # Every omitted source family is checked directly at each new best.
            maximum=-32
            for offset in range(0,len(add),128):
                dot=add[offset:offset+128].astype(np.float32)@x[keep].T.astype(np.float32)
                assert dot.max(initial=-32)<=16;maximum=max(maximum,int(dot.max(initial=-32)))
            dot=add.astype(np.int32)@add.T.astype(np.int32);np.fill_diagonal(dot,-32);assert dot.max(initial=-32)<=16
            np.savez_compressed(a.out/f'coordinates31_{len(f)}.npz',points=2*f,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
            np.savez_compressed(a.out/'best_delta.npz',added=2*add,removed=2*x[remove])
            np.savez_compressed(a.out/'best_heads.npz',points=y,squared_norm=32,inner_product_limit=16)
            best=gain;record.update({'orientation':orientation,'added':len(add),'removed':len(remove),'changed_pair_maximum':maximum,'exact_changed_pairs':'PASS'})
            (a.out/'best.json').write_text(json.dumps(record,indent=2)+'\n')
            print(json.dumps(record),flush=True)
        records.append(record);trial+=1
    report=dict(source=str(a.source),source_sha256=hashlib.sha256(a.source.read_bytes()).hexdigest(),module=str(a.module),module_sha256=hashlib.sha256(a.module.read_bytes()).hexdigest(),trials=len(records),best_gain=best,seconds=time.monotonic()-start,records=records,scope='Exact two-copy union optimization for each sampled E8/D7 orientation. Bounded orientation search, not a global optimum.')
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2),flush=True)


if __name__=='__main__':main()
