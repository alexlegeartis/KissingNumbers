#!/usr/bin/env python3
"""Independently check the repaired lifts from coordinates and seed certificates.

Does not read the construction's pair cache or graph cover. Reconstructs every
retained head's seed index, regenerates all required pair products, and verifies
every forbidden pair is absent. The integer-product absolute-term bound makes
binary64 dot products exact. Group reduction is justified by exact isometries.
"""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'dim73-80'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started=time.time()
    original=json.loads((HERE/'gamma72_coset_certificate.json').read_text())
    for name in ('gamma72_coset_heads.npz','gamma72_coset_certificate.py'):
        assert digest(HERE/name)==original['sha256'][name]
    source=np.load(HERE/'gamma72_coset_heads.npz')
    retained=np.load(HERE/'gamma72_repaired_lifts.npz')
    report=json.loads((HERE/'gamma72_repaired_lifts.json').read_text())
    gram=np.load(SOURCE/'gamma72_gram.npy').astype(np.int64)
    g=np.load(SOURCE/'gamma72_gens.npy')[1].astype(np.int64)
    assert np.array_equal(g@gram@g.T,gram)
    assert np.array_equal(np.linalg.matrix_power(g,12),-np.eye(72,dtype=np.int64))
    powers=[np.linalg.matrix_power(g,k) for k in range(12)]
    inverses=[np.eye(72,dtype=np.int64)]+[-powers[12-k] for k in range(1,12)]
    xs=[];lookups=[];centres=[]
    for h in (12625758,4915087,6896114):
        x=source['head_%d'%h].astype(np.int64)
        chosen=x[np.arange(len(x)),np.argmax(x!=0,axis=1)]>0
        x=x[chosen];xs.append(x)
        lookups.append({tuple(v):j for j,v in enumerate(x)})
        centres.append(x[0]%2)
    requests={};pair_count=0
    for n,entry in report['dimensions'].items():
        tails=retained['n%s_tails'%n].astype(np.int64)
        assignment=retained['n%s_assignment'%n]
        assert tails.shape==(entry['cosets'],int(n)-72)
        assert np.all(np.sum(tails*tails,axis=1)==2)
        tailgram=tails@tails.T
        assert np.all(np.abs(tailgram-np.eye(len(tails),dtype=np.int64)*2)<=1)
        masks=[];parities=[];count=0
        for i,(a,p) in enumerate(assignment):
            x=retained['n%s_heads_%d'%(n,i)].astype(np.int64)
            recovered=x@inverses[p]
            indices=[lookups[a][tuple(v)] for v in recovered]
            assert len(indices)==len(set(indices))
            mask=np.zeros(len(xs[a]),dtype=bool);mask[indices]=True;masks.append(mask)
            assert np.all(np.einsum('ij,ij->i',x@gram,x)==22)
            expected=centres[a]@powers[p]%2
            assert np.all(x%2==expected)
            parities.append(tuple(expected))
            assert len(x)==entry['retained_head_lines'][i]
            count+=4*len(x)
        assert len(set(parities))==len(parities)
        assert 6218175600+count==entry['bound']
        assert entry['gain']==entry['bound']-entry['previous_working_bound']>0
        for i,(a,p) in enumerate(assignment):
            for j in range(i+1,len(assignment)):
                b,q=assignment[j]
                # t_i dot t_j = (5/4) times the integral D_k root dot.
                limit=16-5*abs(int(tailgram[i,j]))
                if a<=b:
                    key=(a,b,int((q-p)%12));left,right=masks[i],masks[j]
                else:
                    key=(b,a,int((p-q)%12));left,right=masks[j],masks[i]
                requests.setdefault(key,[]).append((limit,left,right,n,i,j))
                pair_count+=1
    checked_products=0
    with threadpool_limits(limits=4):
        for number,((a,b,d),checks) in enumerate(sorted(requests.items())):
            left=xs[a]@gram
            right=xs[b]@powers[d]
            term_bound=int(np.max(np.sum(np.abs(left),axis=1)))*int(np.max(np.abs(right)))
            assert term_bound<2**53
            for start in range(0,len(left),400):
                dot=np.abs(left[start:start+400].astype(np.float64)@right.T.astype(np.float64))
                assert np.all(dot==np.rint(dot))
                checked_products+=dot.size
                for limit in set(c[0] for c in checks):
                    rr,cc=np.nonzero(dot>limit)
                    rr+=start
                    for required,lmask,rmask,n,i,j in checks:
                        if required==limit:
                            assert not np.any(lmask[rr]&rmask[cc]),(n,i,j,a,b,d)
            if (number+1)%15==0:print('independent pair families checked',number+1,flush=True)
    # Within each nearest coset, distinct head vectors have distance >=8.
    # Reusing a head at opposite tails has squared distance 10 >=8.
    # Every head stays compatible with all Gamma72 minimal vectors by nearestness.
    report['status']='Verified: explicit heads, exact Gram inequalities, inherited certified nearestness'
    report['verification']={'distinct_coset_pairs':pair_count,
        'independent_product_families':len(requests),
        'exact_representative_products':checked_products,
        'seconds':time.time()-started,
        'within_coset':'lattice minimum 8; antipodal tail distance squared 10',
        'base_shell':'nearestness of the three seed sets and exact lattice automorphisms',
        'coordinate_file_sha256':digest(HERE/'gamma72_repaired_lifts.npz'),
        'seed_certificate_sha256':digest(HERE/'gamma72_coset_certificate.json'),
        'verifier_sha256':digest(Path(__file__))}
    (HERE/'gamma72_repaired_lifts.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['verification'],indent=2),flush=True)


if __name__=='__main__':main()
