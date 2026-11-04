"""Root-sum completion and exact core/tail reflections, with finite-pool replay.

No stored numerical outcome is a proof. A no-gain certificate is an injection
of candidates into distinct incompatible old points, checked independently.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.sparse import csr_matrix,coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.optimize import milp,Bounds,LinearConstraint

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext,digest
from discovery.certification import verify_artifact
from discovery.engines.signed_layers import sign_sqrt3


def require(value,message):
    if not value:raise ValueError(message)


def bad(p,q,threshold):
    """Exact vectorized sign(p+q*sqrt(3)-threshold)>0; small integer inputs."""
    r=p-threshold
    require(max(int(abs(r).max(initial=0)),int(abs(q).max(initial=0)))<1000000,'Squared arithmetic bound')
    return ((r>=0)&(q>=0)&((r!=0)|(q!=0))) | ((r>0)&(q<0)&(r*r>3*q*q)) | ((r<0)&(q>0)&(3*q*q>r*r))


def source_small(seed):
    path=ROOT/f'research/coupled13/seed_{seed}.npz'
    with np.load(path,allow_pickle=False) as z:
        A,B=z['A'].astype(np.int64),z['B'].astype(np.int64)
    artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':A.tolist(),'radical_coordinates':B.tolist()})
    proof=verify_artifact(artifact,RunContext(ROOT))
    require(proof.status=='exact' and len(A)==1154,'Small baseline verification')
    return A,B,path,proof.to_dict()


def generate_small(A,B):
    known={tuple(a)+tuple(b) for a,b in zip(A,B)}
    C,D,origins=[],[],[]
    for at in range(0,len(A),64):
        p=A[at:at+64]@A.T+3*(B[at:at+64]@B.T)
        q=A[at:at+64]@B.T+B[at:at+64]@A.T
        for local,j in zip(*np.where((p==-8)&(q==0))):
            i=at+int(local);j=int(j)
            if j<=i:continue
            c,d=A[i]+A[j],B[i]+B[j]
            key=tuple(c)+tuple(d)
            if key in known:continue
            known.add(key);C.append(c);D.append(d);origins.append([i,j])
    return np.asarray(C,np.int64),np.asarray(D,np.int64),origins


def source_high(dimension):
    package=ROOT/'verifications/improved/dim27-31-joint-refills'
    options=list((package/'data').glob(f'coordinates{dimension}_*.npz'))
    require(len(options)==1,'Unique current high-dimensional source')
    path=options[0]
    spec=importlib.util.spec_from_file_location('_sol_high_geometry',package/'verify_geometry.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    proof=module.verify(path)
    with np.load(path,allow_pickle=False) as z:A=z['points'].astype(np.int64)
    return A,np.zeros_like(A),path,proof


def generate_high(A,seed,reflections,per_map,anchors,per_anchor):
    rng=np.random.default_rng(seed);n,d=A.shape
    require(d>24 and np.all(A%2==0),'Integral sparse reflection precondition')
    known=set(map(tuple,A));C=[];origins=[];maps=[]
    # Four equal coefficients make a Householder denominator 2. Mixing one
    # tail coordinate with three core coordinates leaves the old sextet family.
    for trial in range(reflections):
        ids=sorted(map(int,rng.choice(24,3,replace=False)))+[int(rng.integers(24,d))]
        signs=list(map(int,rng.choice([-1,1],4)))
        normal=np.zeros(d,np.int64);normal[ids]=signs
        matrix=2*np.eye(d,dtype=np.int64)-np.outer(normal,normal)
        require(np.array_equal(matrix@matrix.T,4*np.eye(d,dtype=np.int64)),'Exact mixed isometry')
        indices=rng.choice(n,min(n,per_map),replace=False)
        transformed=A[indices]@matrix
        require(np.all(transformed%2==0),'Integral mixed images')
        transformed//=2
        maps.append({'normal':normal.tolist(),'source_indices':indices.tolist()})
        for source,point in zip(indices,transformed):
            key=tuple(point)
            if key in known:continue
            known.add(key);C.append(point);origins.append(['reflection',trial,int(source)])
    auxiliary=np.flatnonzero(np.any(A[:,24:]!=0,axis=1))
    for anchor in rng.choice(auxiliary,min(anchors,len(auxiliary)),replace=False):
        dots=A@A[anchor]
        possible=np.flatnonzero(dots==-64)
        for j in rng.choice(possible,min(per_anchor,len(possible)),replace=False):
            point=A[anchor]+A[j];key=tuple(point)
            if key in known:continue
            known.add(key);C.append(point);origins.append(['contact',int(anchor),int(j)])
    candidates=np.asarray(C,np.int64)
    require(candidates.ndim==2 and np.all(np.sum(candidates*candidates,axis=1)==128),'Generated high candidate norms')
    return candidates,np.zeros_like(candidates),{'maps':maps,'point_origins':origins,
        'seed':seed,'reflections':reflections,'per_map':per_map,'anchors':anchors,'per_anchor':per_anchor}


def screen(A,B,C,D,norm,cutoff):
    histogram=Counter();keep=[];blockers=[]
    rational_only=not np.any(B) and not np.any(D)
    if rational_only:
        bound=A.shape[1]*int(abs(A).max(initial=0))*int(abs(C).max(initial=0))
        require(bound<2**24,'Exact binary32 candidate dot bound')
        base=A.astype(np.float32)
    for at in range(0,len(C),64):
        c,d=C[at:at+64],D[at:at+64]
        counts=np.zeros(len(c),np.int64);edges=[[] for _ in c]
        for start in range(0,len(A),4096):
            if rational_only:
                p=(c.astype(np.float32)@base[start:start+4096].T).astype(np.int64)
                q=np.zeros_like(p)
            else:
                p=c@A[start:start+4096].T+3*(d@B[start:start+4096].T)
                q=c@B[start:start+4096].T+d@A[start:start+4096].T
            conflicts=bad(p,q,norm//2)
            counts+=conflicts.sum(axis=1)
            for i,row in enumerate(conflicts):
                if counts[i]<=cutoff:edges[i].extend((np.flatnonzero(row)+start).tolist())
        for i,count in enumerate(counts):
            histogram[int(count)]+=1
            if count<=cutoff:
                require(len(edges[i])==count,'Full blocker list')
                keep.append(at+i);blockers.append(edges[i])
    return C[keep],D[keep],blockers,keep,dict(histogram)


def matching(blockers,old_count):
    ptr=np.r_[0,np.cumsum(list(map(len,blockers)))];edges=np.asarray([j for b in blockers for j in b],np.int64)
    graph=csr_matrix((np.ones(len(edges),np.int8),edges,ptr),shape=(len(blockers),old_count))
    owners=maximum_bipartite_matching(graph,perm_type='column')
    return owners,ptr,edges


def optimize(A,B,C,D,blockers,norm,seconds,clique_cuts=False):
    n=len(C);old=sorted(set(j for b in blockers for j in b));oldids={j:n+i for i,j in enumerate(old)}
    rr=[];cc=[];vv=[];upper=[]
    for i,bs in enumerate(blockers):
        for j in bs:
            row=len(upper);rr.extend([row,row]);cc.extend([i,oldids[j]]);vv.extend([1,-1]);upper.append(0)
    candidate_edges=0
    for at in range(0,n,64):
        p=C[at:at+64]@C.T+3*(D[at:at+64]@D.T)
        q=C[at:at+64]@D.T+D[at:at+64]@C.T
        conflicts=bad(p,q,norm//2)
        for local in range(len(p)):
            i=at+local
            for j in np.flatnonzero(conflicts[local,i+1:])+i+1:
                row=len(upper);rr.extend([row,row]);cc.extend([i,int(j)]);vv.extend([1,1]);upper.append(1);candidate_edges+=1
    # A fixed-height binary sign layer has dot = norm - 2*HammingDistance.
    # A radius-one Hamming ball has pair distance <=2, hence is a conflict
    # clique whenever norm>8. These exact structural cuts strengthen the
    # fractional edge model; they do not freeze a layer or old point.
    cuts=0
    if clique_cuts and norm>8 and A.shape[1]-1<=16:
        by_height={};body_width=A.shape[1]-1
        for variable,a,b,constant in (
                [(i,C[i],D[i],0) for i in range(n)] +
                [(oldids[j],A[j],B[j],1) for j in old]):
            if np.any(b) or not np.all(np.abs(a[:-1])==1):continue
            word=sum((int(x)>0)<<bit for bit,x in enumerate(a[:-1]))
            by_height.setdefault(int(a[-1]),{})[word]=(variable,constant)
        seen=set()
        for height,words in by_height.items():
            for center in range(1<<body_width):
                members=[words[w] for w in [center]+[center^(1<<bit) for bit in range(body_width)] if w in words]
                if len(members)<3:continue
                key=tuple(sorted(variable for variable,_ in members))
                if key in seen:continue
                seen.add(key);old_count=sum(constant for _,constant in members)
                require(old_count<=1,'Baseline sign-layer clique')
                row=len(upper)
                for variable,constant in members:
                    rr.append(row);cc.append(variable);vv.append(-1 if constant else 1)
                upper.append(1-old_count);cuts+=1
    matrix=coo_matrix((vv,(rr,cc)),shape=(len(upper),n+len(old))).tocsc()
    r=milp(np.r_[-np.ones(n),np.ones(len(old))],integrality=np.ones(n+len(old)),bounds=Bounds(0,1),
           constraints=LinearConstraint(matrix,-np.inf,upper),options={'time_limit':seconds,'mip_rel_gap':0})
    added=[] if r.x is None else np.flatnonzero(r.x[:n]>.5).tolist()
    removed=[] if r.x is None else [old[j] for j in np.flatnonzero(r.x[n:]>.5)]
    return {'solver_status':int(r.status),'message':str(r.message),'gain':len(added)-len(removed),
            'added':added,'removed':removed,'candidate_conflicts':candidate_edges,
            'exact_binary_ball_clique_cuts':cuts,
            'numerical_gain_upper_bound':-float(r.mip_dual_bound) if getattr(r,'mip_dual_bound',None) is not None else None,
            'scope':'Numerical MILP diagnostics; only full exact output replay can establish an improvement'}


def replay(path):
    """Independent direct scalar arithmetic: no matching/graph solver needed."""
    meta=json.loads((path/'report.json').read_text());source=ROOT/meta['source']
    require(hashlib.sha256(source.read_bytes()).hexdigest()==meta['source_sha256'],'Source hash')
    with np.load(source,allow_pickle=False) as z:
        if 'points' in z.files:A=z['points'].tolist();B=[[0]*len(A[0]) for _ in A]
        else:A=z['A'].tolist();B=z['B'].tolist()
    witness=path/'witness.npz'
    require(hashlib.sha256(witness.read_bytes()).hexdigest()==meta['witness_sha256'],'Witness hash')
    with np.load(witness,allow_pickle=False) as z:C=z['A'].tolist();D=z['B'].tolist();owners=z['owners'].tolist()
    require(len(C)==len(D)==len(owners)==meta['retained_candidates'],'Witness counts')
    require(len(A)==len(B)==meta['old_cardinality'] and bool(A),'Source counts')
    width=len(A[0]);norm=meta['norm_squared']
    require(type(norm) is int and norm>0 and norm%2==0,'Common positive even norm')
    for a,b in zip(A,B):
        require(len(a)==len(b)==width and all(type(x) is int for x in a+b),'Source coordinate shape')
        require(sum(x*y for x,y in zip(a,b))==0 and sum(x*x+3*y*y for x,y in zip(a,b))==norm,'Source exact norm')
    require(len(set(owners))==len(owners) and all(0<=j<len(A) for j in owners),'Injective blocker owners')
    known={tuple(a)+tuple(b) for a,b in zip(A,B)};candidate_keys=set()
    for c,d,j in zip(C,D,owners):
        require(len(c)==len(d)==width and all(type(x) is int for x in c+d),'Candidate coordinate shape')
        key=tuple(c)+tuple(d)
        require(key not in known and key not in candidate_keys,'New distinct candidates')
        candidate_keys.add(key)
        require(sum(x*y for x,y in zip(c,d))==0 and sum(x*x+3*y*y for x,y in zip(c,d))==norm,'Candidate exact norm')
        p=sum(x*y+3*u*v for x,y,u,v in zip(c,A[j],d,B[j]))
        q=sum(x*v+u*y for x,y,u,v in zip(c,A[j],d,B[j]))
        require(sign_sqrt3(p-norm//2,q)>0,'Assigned owner must be incompatible')
    result={'status':'PASS','candidates':len(C),'scope':'Every subset of this explicit candidate pool requires at least as many distinct old deletions; no bound on other candidates'}
    (path/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mode',choices=['13','31','replay'],required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--seed',default='0_2')
    parser.add_argument('--cutoff',type=int,default=16);parser.add_argument('--seconds',type=float,default=30)
    parser.add_argument('--pool-limit',type=int,default=3000)
    parser.add_argument('--clique-cuts',action='store_true')
    args=parser.parse_args()
    if args.mode=='replay':print(json.dumps(replay(args.out)));return
    args.out.mkdir(parents=True,exist_ok=False);start=monotonic()
    if args.mode=='13':
        A,B,source,proof=source_small(args.seed);C,D,provenance=generate_small(A,B);norm=16
    else:
        A,B,source,proof=source_high(31)
        C,D,provenance=generate_high(A,20261104,24,256,32,64);norm=128
    print(json.dumps({'phase':'generated','old':len(A),'new':len(C)}),flush=True)
    generated=len(C)
    C,D,blockers,indices,hist=screen(A,B,C,D,norm,args.cutoff)
    owners,ptr,edges=matching(blockers,len(A))
    saturated=bool(np.all(owners>=0))
    path=args.out/'witness.npz';np.savez_compressed(path,A=C,B=D,owners=owners,blocker_indptr=ptr,blockers=edges)
    report={'mode':args.mode,'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'source_proof':proof,'norm_squared':norm,'old_cardinality':len(A),'generated_candidates':generated,
            'retained_candidates':len(C),'cutoff':args.cutoff,'blocker_histogram':hist,'retained_generated_indices':indices,
            'matching_saturates_candidates':saturated,'witness_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'generation':provenance,'seconds':monotonic()-start,'new_record':False,
            'scope':'Bounded new contact-sum/core-tail candidate grammar; all incumbent points removable. Cutoff exclusions are outside the finite-pool result.'}
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if not saturated and len(C):
        require(0<args.pool_limit<=12000,'Joint pool limit must be 1 through 12000')
        order=sorted(range(len(C)),key=lambda i:(len(blockers[i]),i))[:args.pool_limit]
        subset_C,subset_D=C[order],D[order]
        result=optimize(A,B,subset_C,subset_D,[blockers[i] for i in order],norm,args.seconds,args.clique_cuts)
        result.update(pool_indices=order,full_pool_optimized=len(order)==len(C))
        report['optimization']=result
        if result['gain']>0:
            keep=np.ones(len(A),bool);keep[result['removed']]=False
            AA=np.concatenate([A[keep],subset_C[result['added']]]);BB=np.concatenate([B[keep],subset_D[result['added']]])
            if args.mode=='13':
                artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':AA.tolist(),'radical_coordinates':BB.tolist()})
                report['positive_exact_geometry']=verify_artifact(artifact,RunContext(ROOT)).to_dict()
                (args.out/'positive.json').write_text(json.dumps(artifact.to_dict())+'\n')
            else:
                from importlib.util import spec_from_file_location,module_from_spec
                spec=spec_from_file_location('_sol_positive_high',ROOT/'verifications/improved/dim27-31-joint-refills/verify_geometry.py')
                module=module_from_spec(spec);spec.loader.exec_module(module)
                candidate=args.out/'positive.npz'
                np.savez_compressed(candidate,points=AA.astype(np.int16),squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
                report['positive_exact_geometry']=module.verify(candidate)
    report['seconds']=monotonic()-start
    (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    if saturated:report['independent_replay']=replay(args.out)
    print(json.dumps({k:report[k] for k in ['mode','generated_candidates','retained_candidates','matching_saturates_candidates','seconds','new_record']}),flush=True)


if __name__=='__main__':main()
