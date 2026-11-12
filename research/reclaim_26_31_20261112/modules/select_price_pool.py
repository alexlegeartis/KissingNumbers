"""Materialize complete blocker lists for a selected exact price cutoff."""
import argparse, json, shutil
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('prices',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cutoff',type=int,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    meta=json.loads((a.prices/'report.json').read_text());z=np.load(meta['source']);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:X//=2
    z=np.load(a.prices/'prices.npz');B,T,cost=z['bodies'],z['tails'],z['blocker_counts']
    i,j=np.nonzero(cost<=a.cutoff);C=np.c_[B[i],T[j]];edges=[];ptr=[0]
    for row,expected in zip(C,cost[i,j]):
        dot=np.zeros(len(X),np.int32)
        for k in np.flatnonzero(row):dot+=X[:,k].astype(np.int32)*int(row[k])
        ids=np.flatnonzero(dot>16)
        if len(ids)!=expected:raise ValueError('Price mismatch')
        edges.extend(ids.tolist());ptr.append(len(edges))
    np.savez_compressed(a.out/'pool.npz',candidates=C,blocker_indices=np.asarray(edges,np.int64),blocker_indptr=np.asarray(ptr,np.int64))
    meta.update({'pool_source':str(a.prices),'selected_pool':len(C),'cutoff':a.cutoff,'pool_edges':len(edges),'pool_limit':None,
        'pool_conclusion':'Unresolved selected pool; direct integer blocker lists replayed'})
    (a.out/'report.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps({'candidates':len(C),'edges':len(edges),'old_points':len(np.unique(edges))}),flush=True)
if __name__=='__main__':main()
