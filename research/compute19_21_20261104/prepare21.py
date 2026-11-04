#!/usr/bin/env python3
"""Rebuild Li's 31560-point configuration and the exact unrestricted line-label graph.

Primary source: https://hexagonmath.org/pdf/2610.00248v1, Section 3.
"""
import itertools,json,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;DATA=HERE/'data';ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'compute17_20_20261104'))
from orbits17 import write_graph
from prepare20 import basis_of


def golay_paper():
    B=np.zeros((12,24),dtype=np.int64)
    for i in range(12):
        for j in (0,1,5,6,7,9,11):B[i,i+j]=1
        B[i,23]=B[i,:23].sum()%2
    M=np.arange(4096,dtype=np.int64)
    G=((M[:,None]>>np.arange(12))&1)@B%2
    assert dict(zip(*np.unique(G.sum(1),return_counts=True)))=={0:1,8:759,12:2576,16:759,24:1}
    return G


def eta(word):
    groups=((3,4),(2,11),(2,5,6,10),(1,3,8,14))
    return sum((sum((int(word)>>(p-1))&1 for p in row)%2)<<i for i,row in enumerate(groups))


def prepare():
    DATA.mkdir(exist_ok=True)
    G=golay_paper();powers=np.int64(1)<<np.arange(21)
    words=G[:,:21]@powers
    line_ids=np.flatnonzero(G[:,:21].sum(1)==5)
    lines=words[line_ids]
    assert len(lines)==21
    assert all((int(a)&int(b)).bit_count()==1 for a,b in itertools.combinations(lines,2))
    line_eta=np.asarray([eta(w) for w in lines])
    assert set(line_eta)==set(range(1,16))
    chosen={e:int(np.flatnonzero(line_eta==e)[0]) for e in range(1,16)}
    labels=np.asarray([chosen[eta(w)] if eta(w) else -1 for w in words],dtype=np.int32)
    assert np.sum(labels>=0)==3840
    # Exact verification of the geometric line rule for all 4096*21*21 cases.
    signs=1-2*G[:,:21]
    profiles=np.asarray([[2 if int(L)>>p&1 else 3 for p in range(21)] for L in lines],dtype=np.int64)
    assert np.all((profiles*profiles).sum(1)==164)
    for a in range(21):
        dots=(signs*profiles[a])@profiles.T
        for m in range(4096):
            if m==0:expected=np.ones(21,dtype=bool)
            elif m in line_ids:
                p=int(np.flatnonzero(line_ids==m)[0])
                expected=np.asarray([not(a!=b and a!=p and b!=p) for b in range(21)])
            else:expected=np.zeros(21,dtype=bool)
            assert np.array_equal(dots[m]>82,expected)
    # Every profile fits every fixed base point; check the exact worst signed-octad products.
    blocks=G[(G.sum(1)==8)&(G[:,21:].sum(1)==0),:21]
    assert len(blocks)==210
    for P in profiles:
        for O in blocks:
            values=P[O==1]
            assert (int(values.sum())-2*int(values.min()))**2<=2*164
        assert max(2*(int(P[i])+int(P[j])) for i,j in itertools.combinations(range(21),2))**2<=2*164
    neighbors=np.bitwise_xor.outer(np.arange(4096,dtype=np.int64),line_ids)
    for m in np.flatnonzero(labels>=0):
        for p,n in enumerate(neighbors[m]):
            if labels[n]>=0:assert len({int(labels[m]),int(labels[n]),p})==3
    # Reconstruct the full published baseline as integer rays.
    base=[]
    for i,j in itertools.combinations(range(21),2):
        for si,sj in itertools.product((-2,2),repeat=2):
            v=[0]*21;v[i]=si;v[j]=sj;base.append(v)
    for block in blocks:
        positions=np.flatnonzero(block)
        for bits in range(256):
            if bits.bit_count()%2:
                v=[0]*21
                for j,p in enumerate(positions):v[p]=-1 if bits>>j&1 else 1
                base.append(v)
    base=np.asarray(base,dtype=np.int16)
    selected=np.flatnonzero(labels>=0)
    auxiliary=(profiles[labels[selected]]*signs[selected]).astype(np.int16)
    assert len(base)==27720 and len(auxiliary)==3840
    np.savez_compressed(DATA/'baseline21.npz',base=base,auxiliary=auxiliary,
                        words=words,lines=lines,line_ids=line_ids,profiles=profiles,
                        labels=labels,code=G)
    # Four isomorphic components; a 961-point component would beat the record.
    sb=basis_of(line_ids)
    component=[0]
    for b in sb:component += [m^b for m in component]
    assert len(component)==1024
    component=sorted(component);index={m:i for i,m in enumerate(component)}
    N=np.asarray([[index[m^int(s)] for s in line_ids] for m in component],dtype=np.int32)
    local_labels=labels[component]
    assert np.sum(local_labels>=0)==960
    rows=[]
    for u in range(1024):
        for p in range(21):
            conflict={21*u+q for q in range(21) if q!=p}
            # Its own color must point to a missing word.
            conflict.update(21*int(N[u,p])+q for q in range(21))
            for q in range(21):
                v=int(N[u,q])
                if q!=p:conflict.add(21*v+p)
                # Incoming labels are also forbidden at a selected word.
                conflict.add(21*v+q)
            conflict.discard(21*u+p)
            rows.append(np.asarray(sorted(conflict),dtype=np.uint32))
    edges=write_graph(DATA/'component21.graph',rows)
    seed=np.asarray([21*u+int(p) for u,p in enumerate(local_labels) if p>=0],dtype=np.int32)
    selected_set=set(map(int,seed))
    assert all(selected_set.isdisjoint(map(int,rows[i])) for i in seed)
    np.savez_compressed(DATA/'component21.npz',component=component,neighbors=N,
                        seed=seed,labels=local_labels,line_ids=line_ids)
    with (DATA/'component21.search').open('w') as f:
        f.write(f'{len(rows)} 30600 21\n'+' '.join(['1']*len(rows))+'\n')
        f.write(str(len(seed))+'\n'+' '.join(map(str,seed))+'\n')
    report=dict(dimension=21,baseline=31560,fixed_base=27720,auxiliary=3840,
                component_words=1024,component_selected=960,candidate_labels=len(rows),
                directed_edges=edges,min_degree=min(map(len,rows)),max_degree=max(map(len,rows)),
                exact_rule_checks=4096*21*21,source='https://hexagonmath.org/pdf/2610.00248v1')
    (DATA/'preparation21.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':prepare()
