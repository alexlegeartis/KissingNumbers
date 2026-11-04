"""An alternative 720-bridge coupling from all six K6 factorizations.

For each left edge {u,v}, use the perfect matching shared by the u-th and v-th
factorizations on the right. Every row and column remains a perfect matching.
This is a connected incidence pattern rather than five separate K3,3 blocks.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,itertools,json
from pathlib import Path
import numpy as np
from audit import clebsch_coefficients,leq,FACTORS
ROOT=Path(__file__).resolve().parent


def matchings(vertices):
    if not vertices:yield ();return
    first=vertices[0]
    for j in vertices[1:]:
        for rest in matchings([x for x in vertices if x not in [first,j]]):
            yield (tuple(sorted((first,j))),*rest)


def bridge_supports():
    ms=sorted(set(tuple(sorted(m)) for m in matchings(list(range(6)))))
    edges=sorted(itertools.combinations(range(6),2));factors=[]
    for group in itertools.combinations(range(15),5):
        used=[edge for m in group for edge in ms[m]]
        if len(set(used))==15:factors.append(set(group))
    assert len(ms)==15 and len(factors)==6
    supports=[];matrix=np.zeros((15,15),dtype=np.int64)
    for i,(u,v) in enumerate(edges):
        common=factors[u]&factors[v];assert len(common)==1
        for edge in ms[next(iter(common))]:
            j=edges.index(edge);matrix[i,j]=1;supports.append((u,v,6+edge[0],6+edge[1]))
    assert len(supports)==45 and np.all(matrix.sum(axis=0)==3) and np.all(matrix.sum(axis=1)==3)
    # Each neighborhood is a perfect matching, on both sides.
    for row in [*matrix,*matrix.T]:
        endpoints=[v for i in np.flatnonzero(row) for v in edges[i]]
        assert sorted(endpoints)==list(range(6))
    visited={0};frontier=[0]
    adjacency=[np.flatnonzero(row)+15 for row in matrix]+[np.flatnonzero(row) for row in matrix.T]
    while frontier:
        for v in adjacency[frontier.pop()]:
            if int(v) not in visited:visited.add(int(v));frontier.append(int(v))
    assert len(visited)==30
    return supports,matrix,ms,factors


def coefficients():
    a,b=clebsch_coefficients();supports,_,_,_=bridge_supports();rows=[]
    for support in supports:
        for signs in itertools.product([-1,1],repeat=4):
            row=np.zeros(12,dtype=np.int64);row[list(support)]=6*np.array(signs);rows.append(row)
    a[120:]=rows;b[120:]=0
    return a,b


def graph_moments(a,b):
    pp=a@a.T+2*b@b.T;qq=a@b.T+b@a.T
    valid=leq(pp,qq,72);np.fill_diagonal(valid,True);assert valid.all()
    adj=((pp==72)&(qq==0)).astype(np.int64);np.fill_diagonal(adj,0)
    square=adj@adj
    return dict(contact_pairs=int(adj.sum()//2),trace_A3=int(np.sum(square*adj)),
                trace_A4=int(np.sum(square*square)),
                common_neighbor_histogram={str(k):int(v) for k,v in zip(*np.unique(square,return_counts=True))})


def main():
    supports,matrix,ms,factors=bridge_supports();a,b=coefficients();old_a,old_b=clebsch_coefficients()
    old=graph_moments(old_a,old_b);new=graph_moments(a,b)
    result=dict(exact_kissing_840=True,bridge_points=720,coupling_connected=True,
                standard_coupling_components=5,new_coupling_components=1,
                supports=[list(s) for s in supports],perfect_matchings=ms,
                factorizations=[sorted(f) for f in factors],standard=old,alternative=new,
                contact_graph_moments_different=old!=new)
    (ROOT/'alternative_bridges.json').write_text(json.dumps(result,indent=2)+'\n')
    np.save(ROOT/'runs'/'alternative840.npy',(a+np.sqrt(2)*b)/12)
    print(json.dumps({k:v for k,v in result.items() if k not in ['supports','perfect_matchings','factorizations']},indent=2))


if __name__=='__main__':main()
