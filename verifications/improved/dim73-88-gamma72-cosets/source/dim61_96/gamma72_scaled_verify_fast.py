#!/usr/bin/env python3
"""Verify scaled coset lifts using exact binary32 products.

Every sum of absolute integer terms is checked to be below 2**24. Each
multiplication and every partial sum is therefore represented exactly.
No search pair-profile cache is read.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'dim73-80'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--seeds',type=int,default=4)
    args=ap.parse_args();started=time.time();name='gamma72_scaled_%d'%args.seeds
    report=json.loads((HERE/(name+'.json')).read_text())
    data=np.load(HERE/(name+'.npz'))
    seeddata=np.load(HERE/('gamma72_scaled_seeds_%d.npz'%args.seeds))
    for filename,value in report['source_sha256'].items():assert sha(HERE/filename)==value
    cert=json.loads((HERE/'gamma72_coset_certificate.json').read_text())
    assert sha(HERE/'gamma72_coset_heads.npz')==cert['sha256']['gamma72_coset_heads.npz']
    gram=np.load(SOURCE/'gamma72_gram.npy').astype(np.int64)
    generators=np.load(SOURCE/'gamma72_gens.npy').astype(np.int64)
    order=np.load(HERE/'gamma72_order91.npz');identity=np.eye(72,dtype=np.int64)
    g=identity.copy()
    for j in order['word']:g=g@generators[j]
    assert np.array_equal(g,order['generator']) and np.array_equal(g,seeddata['generator'])
    assert np.array_equal(g@gram@g.T,gram)
    powers=[identity]
    for i in range(91):powers.append(powers[-1]@g)
    assert np.array_equal(powers[-1],identity)
    source=np.load(HERE/'gamma72_coset_heads.npz')
    raw=[]
    for h in (12625758,4915087,6896114):
        x=source['head_%d'%h].astype(np.int64)
        raw.append(x[x[np.arange(len(x)),np.argmax(x!=0,axis=1)]>0])
    xs=[]
    for origin,t in zip(seeddata['origins'],seeddata['transforms']):
        assert np.array_equal(t@gram@t.T,gram)
        xs.append(raw[origin]@t)
    requests={};pairs=0;ranks={}
    for n,entry in report['dimensions'].items():
        tag='n'+n;tails=data[tag+'_tails'].astype(np.int64)
        assignment=data[tag+'_assignment'];q=entry['tail_squared_norm'];m=len(tails)
        assert m==entry['cosets']==len(assignment)
        tg=tails@tails.T
        assert np.all(np.diag(tg)==q)
        assert np.all(2*np.abs(tg-q*np.eye(m,dtype=np.int64))<=q)
        rank=sp.Matrix(tails.T@tails).rank();assert rank==int(n)-72;ranks[n]=rank
        masks=[];parities=[];total=0
        for i,(a,p) in enumerate(assignment):
            indices=data[tag+'_indices_%d'%i]
            assert np.all(indices[1:]>indices[:-1])
            assert np.all(indices>=0) and np.all(indices<len(xs[a]))
            mask=np.zeros(len(xs[a]),dtype=bool);mask[indices]=True;masks.append(mask)
            total+=len(indices)
            parities.append(tuple(xs[a][0]@powers[p]%2))
        assert len(set(parities))==m
        assert total==entry['retained_head_lines']
        assert 6218175600+4*total==entry['bound']
        assert entry['bound']-entry['previous_working_bound']==entry['gain']
        for i,(a,p) in enumerate(assignment):
            for j in range(i+1,m):
                b,r=assignment[j];left,right=masks[i],masks[j]
                if a>b:a0,b0,p0,r0=b,a,r,p;left,right=right,left
                else:a0,b0,p0,r0=a,b,p,r
                d=int((r0-p0)%91)
                if a0==b0 and d>45:d=91-d;left,right=right,left
                bound=(16*q-10*abs(int(tg[i,j])))//q
                requests.setdefault((int(a0),int(b0),d),[]).append((bound,left,right,n,i,j))
                pairs+=1
    products=0
    with threadpool_limits(limits=4):
        for number,((a,b,d),checks) in enumerate(sorted(requests.items())):
            left=xs[a]@gram;right=xs[b]@powers[d]
            assert int(np.max(np.sum(np.abs(left),axis=1)))*int(np.max(np.abs(right)))<2**24
            rt=right.T.astype(np.float32);pieces=[]
            for start in range(0,len(left),400):
                dots=np.abs(left[start:start+400].astype(np.float32)@rt)
                rr,cc=np.nonzero(dots>11)
                vv=dots[rr,cc];assert np.all(vv==np.rint(vv))
                pieces.append(np.column_stack([rr+start,cc,vv]).astype(np.int32))
                products+=dots.size
            bad=np.vstack(pieces)
            for bound,lmask,rmask,n,i,j in checks:
                edge=bad[bad[:,2]>bound]
                assert not np.any(lmask[edge[:,0]]&rmask[edge[:,1]]),(n,i,j,a,b,d)
            if (number+1)%100==0:print('independent families',number+1,'/',len(requests),flush=True)
    report['status']='Verified: exact group actions, tail ranks, and independently regenerated cross-head products'
    report['verification']={'independent_pair_families':len(requests),'distinct_coset_pairs':pairs,
        'exact_representative_products':products,'tail_ranks':ranks,'seconds':time.time()-started,
        'nearestness':'inherited certified nearest sets, preserved by checked integral isometries',
        'internal_compatibility':'minimum Gamma72 squared norm 8; antipodal tails squared distance 10',
        'construction_sha256':sha(HERE/(name+'.npz')),'seed_certificate_sha256':sha(HERE/'gamma72_coset_certificate.json'),
        'verifier_file':Path(__file__).name,'verifier_sha256':sha(Path(__file__)),'gram_sha256':sha(SOURCE/'gamma72_gram.npy')}
    (HERE/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['verification'],indent=2),flush=True)


if __name__=='__main__':main()
