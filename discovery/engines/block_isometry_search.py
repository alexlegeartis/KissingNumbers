"""Jointly synthesize local metric isometries of nonlinear block codes.

This changes the construction, not just its representation: different cells
may use different isometries, and a group assignment may tie selected cells.
Counterexamples generate exact inequalities. The oracle may be incomplete;
only a complete, independent verifier may turn a result into a bound.
"""
from itertools import product
from pathlib import Path
import math, time
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix
from ..core import Artifact, Verification, RunContext, check_params, integer
from .binary_translate_search import _matrix, _weight


def _validate(p,context):
    check_params(p,{'generator','modulus','bit_width','block_map','candidate_offset','references',
                    'reference_offsets','map_menu','groups','minimum_distance','rounds','oracle_steps',
                    'exact_pair_limit','seed','seconds','intersection_probes'})
    q=integer(p.get('modulus'),'modulus',2);width=integer(p.get('bit_width'),'bit_width')
    if q>16 or width>8:raise ValueError('invalid modulus or bit width')
    G=_matrix(p.get('generator'),'generator',high=q);k,n=G.shape;context.limit('max_dimension',width*n);context.limit('max_rank_work',k*n)
    block=p.get('block_map')
    if not isinstance(block,list) or len(block)!=q or any(type(x) is not int or not 0<=x<1<<width for x in block):raise ValueError('invalid block map')
    block=np.asarray(block,dtype=np.int64)
    menu=_matrix(p.get('map_menu'),'map_menu',columns=q,high=1<<width);context.limit('max_codewords',len(menu))
    if len({tuple(row) for row in menu})!=len(menu):raise ValueError('duplicate local maps')
    wt=np.array([int(x).bit_count() for x in range(1<<width)],dtype=np.int64)
    base=wt[block[:,None]^block[None,:]]
    for row in menu:
        if not np.array_equal(wt[row[:,None]^row[None,:]],base):raise ValueError('local map does not preserve all block distances')
    groups=p.get('groups')
    if not isinstance(groups,list) or len(groups)!=n or any(type(x) is not int or x<0 for x in groups):raise ValueError('invalid groups')
    if sorted(set(groups))!=list(range(max(groups)+1)):raise ValueError('groups must be consecutive and nonempty')
    groups=np.asarray(groups,dtype=np.int64);ng=int(groups.max())+1;context.limit('max_codewords',ng*len(menu))
    refs=_matrix(p.get('references'),'references',columns=n,high=1<<width);context.limit('max_codewords',len(refs))
    offsets=_matrix(p.get('reference_offsets',[[0]*n for _ in refs]),'reference_offsets',columns=n,high=q)
    if len(offsets)!=len(refs):raise ValueError('reference offset count mismatch')
    co=_matrix([p.get('candidate_offset',[0]*n)],'candidate_offset',columns=n,high=q)[0]
    target=integer(p.get('minimum_distance'),'minimum_distance')
    if target>n*width:raise ValueError('distance exceeds length')
    return q,width,G,block,menu,groups,refs,offsets,co,target,wt


def search_block_isometries(p,context=None):
    """Search maps and replayable violating word pairs; never certifies a SOTA."""
    context=context or RunContext(Path('.'))
    q,width,G,block,menu,groups,refs,offsets,co,target,wt=_validate(p,context);k,n=G.shape;P=len(menu);ng=int(groups.max())+1
    rounds=integer(p.get('rounds',20),'rounds');context.limit('max_trials',rounds)
    steps=integer(p.get('oracle_steps',500),'oracle_steps');context.limit('max_rank_work',rounds*steps*n*q*len(refs))
    pairlimit=integer(p.get('exact_pair_limit',4096),'exact_pair_limit');context.limit('max_pairs',pairlimit)
    probes=integer(p.get('intersection_probes',8),'intersection_probes');context.limit('max_trials',probes)
    seed=integer(p.get('seed',0),'seed',0);seconds=p.get('seconds',5)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);rng=np.random.default_rng(seed);deadline=time.monotonic()+seconds
    nv=ng*P;rows=[];lo=[];hi=[]
    for g in range(ng):
        row=np.zeros((ng,P));row[g]=1;rows.append(row.ravel());lo.append(1);hi.append(1)
    words=q**k;exact=words*words*len(refs)<=pairlimit
    msgs=np.array(list(product(range(q),repeat=k)),dtype=np.int64) if exact else None
    witnesses=[];candidates=[];status='round_limit';coordinates=np.arange(n)
    for trial in range(rounds):
        remaining=deadline-time.monotonic()
        if remaining<=0:status='time_limit';break
        r=milp(rng.standard_normal(nv),integrality=np.ones(nv),bounds=Bounds(0,1),
          constraints=LinearConstraint(csr_matrix(np.asarray(rows)),lo,hi),options={'time_limit':remaining})
        if r.x is None:status='no_solver_incumbent';break
        choice=np.argmax(r.x.reshape(ng,P),axis=1)
        if np.any(abs(r.x-np.rint(r.x))>1e-6):status='invalid_solver_incumbent';break
        for row,l,u in zip(rows,lo,hi):
            value=sum(row[g*P+choice[g]] for g in range(ng))
            if value<l-1e-7 or value>u+1e-7:raise RuntimeError('invalid MILP incumbent')
        cell_maps=menu[choice[groups]];found=None;tested=0
        if exact:
            aa=(msgs@G+co)%q;images_a=cell_maps[coordinates,aa]
            for ri,ref in enumerate(refs):
                images_b=block[(msgs@G+offsets[ri])%q]
                for i,a in enumerate(msgs):
                    for j,b in enumerate(msgs):
                        tested+=1;dist=_weight(images_a[i]^ref^images_b[j])
                        if dist<target:found=(ri,a.copy(),b.copy(),dist);break
                    if found:break
                if found:break
        else:
            # If the first k columns are triangular with unit diagonal, decode
            # transformed words exactly to detect coincident points at once.
            # This is only a rejection oracle; no misses establish separation.
            inverses=None
            if k<=n and len(set(block.tolist()))==q and all(math.gcd(int(G[i,i]),q)==1 for i in range(k)) and not np.any(np.tril(G[:,:k],-1)):
                inverses=np.full(1<<width,-1,dtype=np.int64);inverses[block]=np.arange(q)
            if inverses is not None:
                for attempt in range(probes):
                    a=np.zeros(k,dtype=np.int64) if attempt==0 else rng.integers(q,size=k);av=(a@G+co)%q
                    candidate=cell_maps[coordinates,av]
                    for ri,ref in enumerate(refs):
                        encoded=inverses[candidate^ref]
                        if np.any(encoded<0):continue
                        b=np.zeros(k,dtype=np.int64)
                        for j in range(k):
                            b[j]=((int(encoded[j])-int(offsets[ri,j])-int(b@G[:,j]))*pow(int(G[j,j]),-1,q))%q
                        if np.array_equal((b@G+offsets[ri])%q,encoded):
                            found=(ri,a.copy(),b.copy(),0);break
                    if found:break
            # Decode an anchored candidate point against a reference layer
            # before the joint search. Each failure is a replayable bad pair.
            anchor_a=np.zeros(k,dtype=np.int64);anchor=cell_maps[coordinates,co]
            ri=int(rng.integers(len(refs)));b=rng.integers(q,size=k);bv=(b@G+offsets[ri])%q
            for step in range(steps//2):
                if found:break
                tested+=1;current=_weight(anchor^refs[ri]^block[bv])
                if current<target:found=(ri,anchor_a.copy(),b.copy(),current);break
                i=int(rng.integers(k));old=int(b[i]);variants=(bv[None,:]+(np.arange(q)[:,None]-old)*G[i])%q
                scores=wt[anchor[None,:]^refs[ri][None,:]^block[variants]].sum(axis=1);temperature=.4+2.6*(1-(step%1000)/1000)
                probs=np.exp(-(scores-scores.min())/temperature);probs/=probs.sum();value=int(rng.choice(q,p=probs));b[i]=value;bv=variants[value]
                if step%1000==999 and len(refs)>1:ri=int(rng.integers(len(refs)));bv=(b@G+offsets[ri])%q
            a=rng.integers(q,size=k);b=rng.integers(q,size=k);ri=int(rng.integers(len(refs)))
            av=(a@G+co)%q;bv=(b@G+offsets[ri])%q
            current=_weight(cell_maps[coordinates,av]^refs[ri]^block[bv])
            for step in range(steps-steps//2):
                if found:break
                tested+=1
                if current<target:found=(ri,a.copy(),b.copy(),current);break
                side=int(rng.integers(2));i=int(rng.integers(k));c=a if side==0 else b;old=int(c[i]);base=av if side==0 else bv
                variants=(base[None,:]+(np.arange(q)[:,None]-old)*G[i])%q
                left=cell_maps[coordinates,variants] if side==0 else cell_maps[coordinates,av][None,:]
                right=block[bv][None,:] if side==0 else block[variants]
                scores=wt[left^refs[ri][None,:]^right].sum(axis=1);temperature=.4+2.6*(1-(step%1000)/1000)
                probs=np.exp(-(scores-scores.min())/temperature);probs/=probs.sum();value=int(rng.choice(q,p=probs));c[i]=value;current=int(scores[value])
                if side==0:av=variants[value]
                else:bv=variants[value]
                if step%1000==999 and len(refs)>1:
                    ri=int(rng.integers(len(refs)));bv=(b@G+offsets[ri])%q;current=_weight(cell_maps[coordinates,av]^refs[ri]^block[bv])
        if found:
            ri,a,b,dist=found;aa=(a@G+co)%q;bb=(b@G+offsets[ri])%q
            if _weight(cell_maps[coordinates,aa]^refs[ri]^block[bb])!=dist or dist>=target:raise RuntimeError('invalid oracle witness')
            cell_cost=wt[menu[:,aa].T^refs[ri,:,None]^block[bb,None]]
            row=np.array([cell_cost[groups==g].sum(axis=0) for g in range(ng)]).ravel()
            rows.append(row);lo.append(target);hi.append(np.inf)
            witnesses.append({'reference':int(ri),'message_a':a.tolist(),'message_b':b.tolist(),'map_choices':choice.tolist(),'distance':int(dist)})
        else:
            candidates.append({'map_choices':choice.tolist(),'exhaustive_cross_check':bool(exact),'pairs_tested':tested})
            row=np.zeros((ng,P));row[np.arange(ng),choice]=1;rows.append(row.ravel());lo.append(-np.inf);hi.append(ng-1)
    return {'parameters':p,'status':status,'oracle':'exhaustive' if exact else 'heuristic','candidates':candidates,'rejection_witnesses':witnesses,
            'scope':'Local maps preserve internal block distances. Parent code distance/cardinality, cross distances and spherical embedding require independent proofs.'}


def block_isometry_search(inputs,params,context):
    if inputs:raise ValueError('block isometry search takes no inputs')
    result=search_block_isometries(params,context)
    return Artifact('search_result',{'backend':'block_isometry_search','format':'block_isometry_search_v1','result':result},{},
                    Verification('unverified','constructive search; independent proof obligations remain',{}))


def verify_artifact(artifact,context=None):
    """Replay all bad-pair witnesses and claimed exhaustive cross checks."""
    context=context or RunContext(Path('.'))
    if artifact.kind!='search_result' or artifact.data.get('backend')!='block_isometry_search' or artifact.data.get('format')!='block_isometry_search_v1':raise ValueError('invalid block search receipt')
    r=artifact.data.get('result');p=r.get('parameters') if isinstance(r,dict) else None
    if not isinstance(p,dict):raise ValueError('missing parameters')
    q,width,G,block,menu,groups,refs,offsets,co,target,wt=_validate(p,context);k,n=G.shape;ng=int(groups.max())+1;idx=np.arange(n)
    w=r.get('rejection_witnesses');c=r.get('candidates')
    if not isinstance(w,list) or not isinstance(c,list):raise ValueError('invalid evidence lists')
    context.limit('max_trials',len(w)+len(c))
    def vec(v,length,maximum):
        if not isinstance(v,list) or len(v)!=length or any(type(x) is not int or not 0<=x<maximum for x in v):raise ValueError('invalid evidence vector')
        return np.asarray(v,dtype=np.int64)
    for x in w:
        if not isinstance(x,dict) or type(x.get('reference')) is not int or not 0<=x['reference']<len(refs):raise ValueError('invalid witness reference')
        a=vec(x.get('message_a'),k,q);b=vec(x.get('message_b'),k,q);choices=vec(x.get('map_choices'),ng,len(menu));ri=x['reference']
        dist=_weight(menu[choices[groups],(a@G+co)%q]^refs[ri]^block[(b@G+offsets[ri])%q])
        if type(x.get('distance')) is not int or x['distance']!=dist or dist>=target:raise ValueError('invalid rejection witness')
    exact=[x for x in c if isinstance(x,dict) and x.get('exhaustive_cross_check') is True]
    context.limit('max_pairs',(q**k)**2*len(refs)*len(exact))
    if exact:
        msgs=np.array(list(product(range(q),repeat=k)),dtype=np.int64);aa=(msgs@G+co)%q
        for x in exact:
            choices=vec(x.get('map_choices'),ng,len(menu));images_a=menu[choices[groups],aa]
            for ri,ref in enumerate(refs):
                images_b=block[(msgs@G+offsets[ri])%q]
                for a in images_a:
                    for b in images_b:
                        if _weight(a^ref^b)<target:raise ValueError('false exhaustive cross-distance claim')
    return Verification('unverified','replayed isometry checks, rejection witnesses and claimed exhaustive cross checks',
                        {'facts':{},'witnesses_checked':len(w),'cross_checks_replayed':len(exact)})

OPERATORS={'search.block_isometries':block_isometry_search}
INPUT_KINDS={'search.block_isometries':()}
OUTPUT_KINDS={'search.block_isometries':('search_result',)}
