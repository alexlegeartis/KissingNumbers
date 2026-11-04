"""Build the exact conflict graph of 20736 eight-point quaternionic sign orbits."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np
from numba import njit,prange,set_num_threads
from exact_cells import exact_frames,scalar_table
from search import ROOT


@njit(cache=True,parallel=True)
def conflict_graph(nodes,table):
    n=len(nodes);bad=np.zeros((n,n),dtype=np.uint8)
    for i in prange(n):
        h,f,g,k,a,b,c=nodes[i]
        for j in range(i):
            hh,ff,gg,kk,aa,bb,cc=nodes[j]
            ra=0;rb=0
            for block in range(3):
                if block==0:qa=table[f,ff,a,aa,0];qb=table[f,ff,a,aa,1]
                elif block==1:qa=table[g,gg,b,bb,0];qb=table[g,gg,b,bb,1]
                else:qa=table[k,kk,c,cc,0];qb=table[k,kk,c,cc,1]
                if h==hh:
                    wa=2 if h==block else 1;wb=0
                elif h==block or hh==block:wa=0;wb=1
                else:wa=1;wb=0
                ra+=wa*qa+2*wb*qb;rb+=wa*qb+wb*qa
            if ra>4 or (ra==4 and rb>0) or (ra<4 and 2*rb*rb>(4-ra)*(4-ra)):
                bad[i,j]=1;bad[j,i]=1
    return bad


def main():
    p=argparse.ArgumentParser();p.add_argument('--threads',type=int,default=96)
    a=p.parse_args();set_num_threads(a.threads);start=time.monotonic()
    out=ROOT/'runs'/'orbit_pool';out.mkdir(parents=True,exist_ok=True)
    c,d=exact_frames();table=scalar_table(c,d)
    nodes=[]
    for h in range(3):
        for f,g,k in itertools.product(range(6),repeat=3):
            if [f,g,k][h]<3:continue
            for aa,bb,cc in itertools.product(range(4),repeat=3):nodes.append([h,f,g,k,aa,bb,cc])
    nodes=np.array(nodes,dtype=np.int64);lookup={tuple(row):i for i,row in enumerate(nodes)}
    bad=conflict_graph(nodes,table)
    report=json.loads((ROOT/'runs'/'quaternion_cells'/'report.json').read_text())
    cells=report['cells'];seeds=[]
    for code in report['codes']:
        selected=[]
        for cellid in code['cells']:
            h,f,g,k,delta=cells[cellid]
            selected.extend(lookup[(h,f,g,k,aa,bb,aa^bb^delta)] for aa,bb in itertools.product(range(4),repeat=2))
        assert len(selected)==96 and not bad[np.ix_(selected,selected)].any()
        seeds.append(selected)
    # Exact coordinates (A+B sqrt(2))/4 for one representative of each orbit.
    aa=np.zeros((len(nodes),12),dtype=np.int64);bb=np.zeros_like(aa)
    for i,(h,f,g,k,a,b,v) in enumerate(nodes):
        for j,(frame,label) in enumerate(zip([f,g,k],[a,b,v])):
            if j==h:aa[i,4*j:4*j+4]=2*d[frame-3,label]
            elif frame<3:aa[i,4*j:4*j+4]=c[frame,label]
            else:bb[i,4*j:4*j+4]=d[frame-3,label]
    assert np.all(np.sum(aa*aa+2*bb*bb,axis=1)==16)
    reps=(aa+np.sqrt(2)*bb)/4
    np.save(out/'nodes.npy',nodes);np.savez(out/'exact_representatives.npz',A=aa,B=bb)
    bad.tofile(out/'conflicts.bin')
    np.array([len(nodes),len(seeds)],dtype=np.uint32).tofile(out/'header.bin')
    np.array(seeds,dtype=np.uint16).tofile(out/'seeds.bin')
    reps.astype(np.float64).tofile(out/'representatives.bin')
    result=dict(nodes=len(nodes),conflicting_pairs=int(bad.sum()//2),whole_cell_seeds=len(seeds),
                graph_sha256=hashlib.sha256(bad.tobytes()).hexdigest(),
                arithmetic='Exact comparison of (a+b sqrt(2))/8 to 1/2',
                target_orbits=97,target_points=848,axial_points=72,
                elapsed_seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
