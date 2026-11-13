"""Bounded exact-cover search for a 15D module containing an E8 equator.

This explores a structural change of the Leech--Sloane 2564-point recipe.
No claim follows unless the explicit coordinates pass complete pair checking.
"""
import argparse
from itertools import combinations,product
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix


def model():
    triples=[]
    for groups in combinations(range(7),3):
        for signs in product([0,1],repeat=3):
            triples.append(tuple(2*g+s for g,s in zip(groups,signs)))
    lookup={t:i for i,t in enumerate(triples)}
    blocks=[];rows=[];cols=[]
    for groups in combinations(range(7),4):
        for signs in product([0,1],repeat=4):
            b=tuple(2*g+s for g,s in zip(groups,signs));j=len(blocks);blocks.append(b)
            for t in combinations(b,3):rows.append(lookup[t]);cols.append(j)
    a=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(triples),len(blocks))).tocsc()
    return blocks,a


def parity_basis(words):
    origin=int(words[0]);basis={}
    for value in words:
        w=int(value)^origin
        while w:
            k=w.bit_length()-1
            if k not in basis:basis[k]=w;break
            w^=basis[k]
    assert len(words)==2**len(basis)
    return origin,list(basis.values())


def aligned_pairs(source):
    y=np.load(source)
    words=np.sum((y[y[:,-1]==-1,:14]<0).astype(np.int64)*(1<<np.arange(14)),axis=1)
    _,basis=parity_basis(words)
    dual=[w for w in range(1,1<<14) if all((w&b).bit_count()%2==0 for b in basis)]
    dual_basis={}
    for w0 in dual:
        w=w0
        while w:
            k=w.bit_length()-1
            if k not in dual_basis:dual_basis[k]=w;break
            w^=dual_basis[k]
    checks=list(dual_basis.values());assert len(checks)==5
    labels=[sum(((w>>j)&1)<<k for k,w in enumerate(checks)) for j in range(14)]
    labelset=set(labels)
    offsets=[d for d in range(1,32) if {v^d for v in labels}==labelset]
    assert len(offsets)==1
    d=offsets[0];pairs=[];seen=set()
    for i,v in enumerate(labels):
        if i in seen:continue
        j=labels.index(v^d);pairs.append((i,j));seen.update([i,j])
    assert len(pairs)==7
    return pairs


def build_coordinates(source,transversal_blocks):
    y=np.load(source).astype(np.int16)
    pairs=aligned_pairs(source)
    # Preserve the original two 512-word glue codes, in the paired order.
    order=[i for pair in pairs for i in pair]
    glue=y[np.abs(y[:,-1])==1].copy()
    glue[:,:14]=glue[:,order]
    rows=[]
    for i in range(14):
        for sign in [-4,4]:
            r=np.zeros(15,np.int16);r[i]=sign;rows.append(r)
    quads=[tuple(sorted((2*i,2*i+1,2*j,2*j+1))) for i,j in combinations(range(7),2)]
    quads += list(transversal_blocks)
    assert len({t for b in quads for t in combinations(b,3)})==4*len(quads)
    for support in quads:
        for signs in product([-2,2],repeat=4):
            r=np.zeros(15,np.int16);r[list(support)]=signs;rows.append(r)
    rows.extend(glue)
    for i in range(7):
        for a,b,c in product([-2,2],repeat=3):
            r=np.zeros(15,np.int16);r[2*i]=a;r[2*i+1]=b;r[14]=c;rows.append(r)
    yy=np.array(rows,np.int16)
    z=np.zeros_like(yy)
    # First seven coordinates are pair sums; eighth is the A1 coordinate.
    # Last seven coordinates are pair differences.
    for i in range(7):z[:,i]=yy[:,2*i]+yy[:,2*i+1];z[:,8+i]=yy[:,2*i]-yy[:,2*i+1]
    z[:,7]=2*yy[:,14]
    assert np.all(np.sum(z.astype(np.int32)**2,axis=1)==32)
    assert len(z)==len(np.unique(z,axis=0))
    gram=z.astype(np.int32)@z.T.astype(np.int32);np.fill_diagonal(gram,-32)
    assert int(gram.max())<=16
    equator=~np.any(z[:,8:],axis=1)
    return z,pairs,int(np.sum(equator))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--seconds',type=float,default=25)
    ap.add_argument('--exact',action='store_true')
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    blocks,m=model();lower=np.zeros(len(blocks));upper=np.ones(len(blocks));lower[0]=1
    r=milp(-np.ones(len(blocks)),integrality=np.ones(len(blocks)),bounds=Bounds(lower,upper),
           constraints=LinearConstraint(m,np.ones(m.shape[0]) if a.exact else np.zeros(m.shape[0]),np.ones(m.shape[0])),
           options={'time_limit':a.seconds,'presolve':False})
    report={'status':int(r.status),'message':r.message,'seconds':time.monotonic()-started,
            'exact_cover_requested':a.exact,'candidate_blocks':len(blocks),'triples':m.shape[0],
            'numerical_dual_bound':getattr(r,'mip_dual_bound',None)}
    if r.x is not None:
        selected=np.flatnonzero(r.x>0.5);assert np.all(m[:,selected].sum(axis=1)<=1)
        chosen=[blocks[i] for i in selected];report['transversal_blocks']=len(chosen)
        report['blocks']=chosen
        z,pairs,equator=build_coordinates(a.source,chosen)
        report.update({'points':len(z),'aligned_pairs':pairs,'equator_points':equator,
                       'tail_norm_histogram':{str(int(k)):int(v) for k,v in zip(*np.unique(np.sum(z[:,8:].astype(np.int32)**2,axis=1),return_counts=True))},
                       'full_integer_geometry':'PASS'})
        np.savez_compressed(a.out/'module.npz',points=z,squared_norm=32,inner_product_limit=16)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
