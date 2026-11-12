"""Exact finite two-copy exchanges with Barnes--Wall normal modules.

The 15-dimensional coordinate section of BW16 has 2340 norm32 points:
420 roots (±4,±4) and 15 affine-hyperplane supports, each with 128 even signs.
Embed its eight-coordinate E8 equator into a Golay octad.  Reorient the E8
factor and the seven tail coordinates, and select a compatible union by exact
bipartite matching.  Baseline validity must be independently replayed before
calling a positive output a lower bound.
"""
import argparse, hashlib, json, time
from pathlib import Path
from itertools import combinations, product
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

NORMAL=np.asarray([0,2,4,5,6,10,11,23])

def bw15():
    rows=[]
    for i,j in combinations(range(15),2):
        for a,b in product((-4,4),repeat=2):
            x=np.zeros(15,np.int16);x[i]=a;x[j]=b;rows.append(x)
    supports=[]
    for a in range(1,16):
        value=(a&15).bit_count()%2 ^ 1
        support=[i for i in range(16) if ((a&i).bit_count()%2)==value]
        assert len(support)==8 and 15 not in support
        supports.append(support)
        for signs in product((-2,2),repeat=8):
            if sum(s<0 for s in signs)%2:continue
            x=np.zeros(15,np.int16);x[support]=signs;rows.append(x)
    X=np.asarray(rows,np.int16)
    assert len(X)==2340 and len(np.unique(X,axis=0))==2340
    assert np.all(np.sum(X.astype(np.int64)**2,axis=1)==32)
    D=X.astype(np.int32)@X.astype(np.int32).T;np.fill_diagonal(D,-32)
    assert D.max()==16
    q=np.sum(X[:,8:].astype(np.int64)**2,axis=1)
    assert dict(zip(*np.unique(q,return_counts=True)))=={0:240,16:2016,32:84}
    return X[q==16]

def orientation(Y,trial,rng):
    Z=Y.copy();spec={}
    if trial:
        perm=rng.permutation(8);signs=rng.choice([-1,1],8);signs[-1]=np.prod(signs[:-1])
        Z[:,:8]=Z[:,perm]*signs;spec['core_permutation']=perm.tolist();spec['core_signs']=signs.tolist()
        if trial%3==1:
            p=rng.permutation(8);blocks=[p[:4],p[4:]]
            for b in blocks:
                assert np.all(Z[:,b].sum(axis=1)%2==0)
                Z[:,b]=Z[:,b].sum(axis=1)[:,None]//2-Z[:,b]
            Z[:,blocks[0]]*=-1;spec['core_tetrads']=[b.tolist() for b in blocks]
        roots=[]
        for _ in range(trial%5):
            s=rng.choice([-1,1],8);s[-1]=np.prod(s[:-1]);dots=Z[:,:8]@s
            assert np.all(dots%4==0)
            Z[:,:8]-=(dots//4)[:,None]*s;roots.append(s.tolist())
        spec['core_root_reflections']=roots
        p=rng.permutation(7);s=rng.choice([-1,1],7);Z[:,8:]=Z[:,8+p]*s
        spec['tail_permutation']=p.tolist();spec['tail_signs']=s.tolist()
    assert np.all(np.sum(Z[:,:8].astype(np.int64)**2,axis=1)==16)
    assert np.all(np.sum(Z[:,8:].astype(np.int64)**2,axis=1)==16)
    return Z,spec

def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--trials',type=int,default=64);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();z=np.load(args.source,allow_pickle=False);X=z['points'].astype(np.int16)
    if int(z['squared_norm'])==128:
        assert not np.any(X%2);X//=2
    assert X.shape[1]==31 and np.all(np.sum(X.astype(np.int64)**2,axis=1)==32)
    q=np.sum(X[:,24:].astype(np.int64)**2,axis=1);h=np.sum(X[:,NORMAL].astype(np.int64)**2,axis=1)
    assert set(q.tolist())=={0,8,16,32}
    proj=X[:,np.r_[NORMAL,np.arange(24,31)]]
    relevant=((q==8)&(h>8))|((q==16)&(h>0));ids=np.flatnonzero(relevant)
    module_mask=~np.any(X[:,[j for j in range(24) if j not in NORMAL]],axis=1)
    external_columns=~module_mask[ids]
    safe_core=np.unique(proj[q==0,:8],axis=0).astype(np.float32)
    safe_mid=np.unique(proj[(q==8)&~relevant,:8],axis=0).astype(np.float32)
    mid_tails=np.unique(proj[(q==8)&~relevant,8:],axis=0).astype(np.float32)
    pure_tails=np.unique(proj[q==32,8:],axis=0).astype(np.float32)
    old=proj[ids];oldset={r.tobytes():i for i,r in enumerate(proj)}
    base=bw15();rng=np.random.default_rng(20261112);records=[];bestgain=0
    print(json.dumps({'baseline':len(X),'normal_module':int(np.sum(np.sum(X[:,[j for j in range(24) if j not in NORMAL]]**2,axis=1)==0)),
        'relevant_old':len(ids),'midpoint_exceptions':int(np.sum(relevant&(q==8))),
        'head_exceptions_including_module':int(np.sum(relevant&(q==16))),'new_module_heads':len(base)}),flush=True)
    for trial in range(args.trials):
        Y,spec=orientation(base,trial,rng)
        yf=Y.astype(np.float32);bc=yf[:,:8];bt=yf[:,8:]
        # These finite exact-integer guards cover every excluded old point.
        assert np.max(bc@safe_core.T)<=16
        assert np.max(bc@safe_mid.T+0,axis=1).max()<=8
        assert np.max(bt@mid_tails.T,axis=1).max()<=8
        assert np.max(bt@pure_tails.T)<=16
        new=np.asarray([r for r in Y if r.tobytes() not in oldset],np.int16).reshape((-1,15))
        cross=new.astype(np.float32)@old.astype(np.float32).T>16
        graph=csr_matrix(cross)
        match=maximum_bipartite_matching(graph,perm_type='column')
        external_count=int(np.sum(np.any(cross[:,external_columns],axis=0)))
        gain=len(new)-int(np.sum(match>=0));record={'trial':trial,'gain':gain,'new_points':len(new),'cross_conflicts':graph.nnz,
            'external_blocker_count':external_count,'whole_module_replacement_gain':58-external_count,'orientation':spec}
        if gain>0:
            rev=np.full(len(old),-1,np.int64)
            for i,j in enumerate(match):
                if j>=0:rev[j]=i
            left=set(map(int,np.flatnonzero(match<0)));right=set();queue=list(left)
            for i in queue:
                for j0 in graph.indices[graph.indptr[i]:graph.indptr[i+1]]:
                    j=int(j0)
                    if match[i]==j or j in right:continue
                    right.add(j)
                    if rev[j]>=0 and int(rev[j]) not in left:left.add(int(rev[j]));queue.append(int(rev[j]))
            remove=ids[sorted(right)];add15=new[sorted(left)];add=np.zeros((len(add15),31),np.int16);add[:,NORMAL]=add15[:,:8];add[:,24:]=add15[:,8:]
            keep=np.ones(len(X),bool);keep[remove]=False;F=np.r_[X[keep],add]
            assert len(F)==len(X)+gain and len(np.unique(F,axis=0))==len(F)
            # Independent direct integer changed-pair verification.
            for i,row in enumerate(add):
                products=np.zeros(len(F),np.int32)
                for j in np.flatnonzero(row):products+=F[:,j].astype(np.int32)*int(row[j])
                products[len(X)-len(remove)+i]=-32
                assert products.max(initial=-32)<=16
            record.update({'added':len(add),'removed':len(remove),'exact_changed_pairs':'PASS'})
            np.savez_compressed(args.out/f'proposal31_{len(F)}_trial{trial}.npz',points=2*F,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
            np.savez_compressed(args.out/f'exchange_trial{trial}.npz',removed=X[remove],added=add)
            bestgain=max(bestgain,gain)
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k!='orientation'}),flush=True)
        report={'source':str(args.source.resolve()),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'source_count':len(X),
            'trials_requested':args.trials,'trials_completed':len(records),'best_gain':bestgain,'seconds':time.monotonic()-start,
            'scope':'Exact finite two-copy unions. New BW15 module and changed pairs are exact; source full geometry remains independently required.',
            'relevant_old':len(ids),'normal_octad':NORMAL.tolist(),'records':records}
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'best_gain':bestgain,'seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__':main()
