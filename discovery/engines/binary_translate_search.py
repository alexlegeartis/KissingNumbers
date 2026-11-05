"""Synthesize binary translates of nonlinear block images using real counterexamples.

The underlying ring code is represented by a generator, but XOR in the binary
image is *not* replaced by ring addition. A failed proposal yields a pair of
ring messages and a reference layer, which produces a valid inequality for all
future proposals. Numerical MILP only proposes translates. A bounded local
oracle can reject a proposal but cannot certify its acceptance.
"""
from itertools import product
from pathlib import Path
import math
import time
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix
from ..core import Artifact, Verification, RunContext, integer, check_params


def _matrix(v, name, columns=None, low=0, high=None):
    if not isinstance(v, list) or not v or any(not isinstance(row,list) or not row for row in v):
        raise ValueError(f'{name} must be a nonempty integer matrix')
    width=len(v[0])
    if (columns is not None and width!=columns) or any(len(row)!=width for row in v):
        raise ValueError(f'{name} has inconsistent width')
    if any(type(x) is not int or x<low or (high is not None and x>=high) for row in v for x in row):
        raise ValueError(f'{name} contains invalid integers')
    return np.asarray(v,dtype=np.int64)


def _weight(x):
    return sum(int(v).bit_count() for v in x)


def _difference(a,b,G,q,block_map,offset_a,offset_b):
    return block_map[(a@G+offset_a)%q] ^ block_map[(b@G+offset_b)%q]


def search_binary_translates(params, context=None):
    """Produce exact counterexamples and candidate translates, not a record.

    ``alphabet`` lists permitted binary block values independently per cell;
    ``references`` are previously chosen binary translates. The searchable
    base is {block_map[a G] : a in Z_q^k}, without assuming an additive Gray
    map. Output ``exhaustive_cross_check`` only concerns cross-distances to
    references; the base code's internal distance/cardinality remain separate
    proof obligations.
    """
    context=context or RunContext(Path('.'))
    allowed={'generator','modulus','block_map','bit_width','alphabet','references',
             'minimum_distance','rounds','oracle_steps','exact_pair_limit','seed','seconds',
             'candidate_offset','reference_offsets'}
    check_params(params,allowed)
    q=integer(params.get('modulus'),'modulus',2)
    if q>16:raise ValueError('modulus must be at most 16')
    width=integer(params.get('bit_width'),'bit_width')
    if width>8:raise ValueError('bit_width must be at most 8')
    block=params.get('block_map')
    if not isinstance(block,list) or len(block)!=q or any(type(x) is not int or not 0<=x<2**width for x in block):
        raise ValueError('block_map must provide one valid binary block per ring symbol')
    block=np.asarray(block,dtype=np.int64)
    G=_matrix(params.get('generator'),'generator',high=q);k,n=G.shape
    context.limit('max_dimension',n*width);context.limit('max_rank_work',k*n)
    refs=_matrix(params.get('references'),'references',columns=n,high=2**width);context.limit('max_codewords',len(refs))
    ring_offsets=_matrix(params.get('reference_offsets',[[0]*n for _ in refs]),'reference_offsets',columns=n,high=q)
    if len(ring_offsets)!=len(refs):raise ValueError('reference_offsets must match references')
    candidate_offset=_matrix([params.get('candidate_offset',[0]*n)],'candidate_offset',columns=n,high=q)[0]
    alphabet=params.get('alphabet')
    if not isinstance(alphabet,list) or len(alphabet)!=n or any(not isinstance(x,list) or not x for x in alphabet):
        raise ValueError('alphabet needs a nonempty list per cell')
    if any(any(type(v) is not int or not 0<=v<2**width for v in a) or len(set(a))!=len(a) for a in alphabet):
        raise ValueError('alphabet blocks must be distinct valid integers')
    target=integer(params.get('minimum_distance'),'minimum_distance')
    if target>n*width:raise ValueError('distance exceeds binary length')
    rounds=integer(params.get('rounds',20),'rounds');context.limit('max_trials',rounds)
    steps=integer(params.get('oracle_steps',500),'oracle_steps')
    context.limit('max_rank_work',rounds*steps*n*q*len(refs))
    pair_limit=integer(params.get('exact_pair_limit',4096),'exact_pair_limit');context.limit('max_pairs',pair_limit)
    seed=integer(params.get('seed',0),'seed',0)
    seconds=params.get('seconds',5)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds)
    rng=np.random.default_rng(seed);deadline=time.monotonic()+seconds
    offsets=np.cumsum([0]+[len(a) for a in alphabet]);nv=int(offsets[-1]);context.limit('max_codewords',nv)
    wt=np.array([i.bit_count() for i in range(1<<width)],dtype=np.int64)
    rows=[];lower=[];upper=[]
    for j in range(n):
        row=np.zeros(nv);row[offsets[j]:offsets[j+1]]=1;rows.append(row);lower.append(1);upper.append(1)
    words=q**k
    exact=words*words*len(refs)<=pair_limit
    all_messages=None
    if exact:all_messages=np.array(list(product(range(q),repeat=k)),dtype=np.int64)
    witnesses=[];candidates=[];outcome='round_limit'
    for trial in range(rounds):
        remaining=deadline-time.monotonic()
        if remaining<=0:outcome='time_limit';break
        result=milp(rng.standard_normal(nv),integrality=np.ones(nv),bounds=Bounds(np.zeros(nv),np.ones(nv)),
          constraints=LinearConstraint(csr_matrix(np.asarray(rows)),lower,upper),options={'time_limit':remaining})
        if result.x is None:
            outcome='no_solver_incumbent';break
        choice=[int(np.argmax(result.x[offsets[j]:offsets[j+1]])) for j in range(n)]
        proposal=np.array([alphabet[j][choice[j]] for j in range(n)],dtype=np.int64)
        # An incumbent is accepted only after checking its discrete constraints.
        if np.any(np.abs(result.x-np.rint(result.x))>1e-6) or any(sum(result.x[offsets[j]:offsets[j+1]])<.999 for j in range(n)):
            outcome='invalid_solver_incumbent';break
        for row,lo,hi in zip(rows,lower,upper):
            z=sum(row[offsets[j]+choice[j]] for j in range(n))
            if z<lo-1e-7 or z>hi+1e-7:raise RuntimeError('MILP incumbent violates proposal constraints')
        found=None;tested=0
        if exact:
            # Enumerate messages independently, including a=b and duplicates.
            images_a=block[(all_messages@G+candidate_offset)%q]
            for ri,ref in enumerate(refs):
                images_b=block[(all_messages@G+ring_offsets[ri])%q]
                for i,a in enumerate(all_messages):
                    for j,b in enumerate(all_messages):
                        tested+=1;d=images_a[i]^images_b[j];dist=_weight(proposal^ref^d)
                        if dist<target:found=(ri,a.copy(),b.copy(),d,dist);break
                    if found:break
                if found:break
        else:
            # Local search for an explicit bad pair. Failure to find one is not proof.
            a=rng.integers(q,size=k);b=rng.integers(q,size=k)
            ri=int(rng.integers(len(refs)));av=(a@G+candidate_offset)%q;bv=(b@G+ring_offsets[ri])%q
            current=_weight(proposal^refs[ri]^block[av]^block[bv])
            for step in range(steps):
                tested+=1
                if current<target:
                    d=block[av]^block[bv];found=(ri,a.copy(),b.copy(),d,current);break
                side=int(rng.integers(2));i=int(rng.integers(k));c=a if side==0 else b;old=int(c[i])
                # Explore all values of one message coordinate at once. This
                # is a proposal oracle only, not an exact distance algorithm.
                values=np.arange(q,dtype=np.int64);base=av if side==0 else bv
                changes=(base[None,:]+(values[:,None]-old)*G[i])%q
                xor=(proposal^refs[ri])[None,:]^(block[changes] if side==0 else block[av][None,:])^(block[bv][None,:] if side==0 else block[changes])
                scores=wt[xor].sum(axis=1);temperature=0.4+2.6*(1-(step%1000)/1000)
                probabilities=np.exp(-(scores-scores.min())/temperature);probabilities=probabilities/probabilities.sum()
                value=int(rng.choice(q,p=probabilities));c[i]=value;current=int(scores[value])
                if side==0:av=changes[value]
                else:bv=changes[value]
                if step%1000==999 and len(refs)>1:
                    ri=int(rng.integers(len(refs)));bv=(b@G+ring_offsets[ri])%q
                    current=_weight(proposal^refs[ri]^block[av]^block[bv])
        if found:
            ri,a,b,d,dist=found
            assert np.array_equal(d,_difference(a,b,G,q,block,candidate_offset,ring_offsets[ri])) and _weight(proposal^refs[ri]^d)==dist and dist<target
            row=np.concatenate([wt[np.asarray(alphabet[j])^refs[ri,j]^d[j]] for j in range(n)]).astype(float)
            rows.append(row);lower.append(target);upper.append(np.inf)
            witnesses.append({'reference':int(ri),'message_a':a.tolist(),'message_b':b.tolist(),'rejected_translate':proposal.tolist(),'distance':int(dist)})
        else:
            candidates.append({'translate':proposal.tolist(),'exhaustive_cross_check':bool(exact),'pairs_tested':tested})
            # Distinct proposal next time; cross-certification does not certify internal distance.
            row=np.zeros(nv)
            for j in range(n):row[offsets[j]+choice[j]]=1
            rows.append(row);lower.append(-np.inf);upper.append(n-1)
    return {'parameters':params,'status':outcome,'candidates':candidates,'rejection_witnesses':witnesses,
            'searches':len(witnesses)+len(candidates),'oracle':'exhaustive' if exact else 'heuristic',
            'scope':'Cross-layer tests only. Base-code internal distance, cardinality and spherical embedding require independent proof.'}


def binary_translate_search(inputs,params,context):
    if inputs:raise ValueError('binary translate search takes no inputs')
    data=search_binary_translates(params,context)
    return Artifact('search_result',{'backend':'binary_translate_search','format':'binary_translate_search_v1','result':data},{},
                    Verification('unverified','constructive search; parent and spherical proof obligations remain',{}))


def verify_artifact(artifact, context=None):
    """Replay claimed negative witnesses and any exhaustive cross checks only."""
    context=context or RunContext(Path('.'))
    if artifact.kind!='search_result' or artifact.data.get('backend')!='binary_translate_search' or artifact.data.get('format')!='binary_translate_search_v1':
        raise ValueError('not a binary translate search receipt')
    result=artifact.data.get('result');p=result.get('parameters') if isinstance(result,dict) else None
    if not isinstance(p,dict):raise ValueError('missing search parameters')
    q=integer(p.get('modulus'),'modulus',2);width=integer(p.get('bit_width'),'bit_width')
    if q>16 or width>8:raise ValueError('invalid modulus or bit width')
    G=_matrix(p.get('generator'),'generator',high=q);k,n=G.shape;context.limit('max_dimension',n*width);context.limit('max_rank_work',k*n)
    refs=_matrix(p.get('references'),'references',columns=n,high=1<<width)
    ring_offsets=_matrix(p.get('reference_offsets',[[0]*n for _ in refs]),'reference_offsets',columns=n,high=q)
    if len(ring_offsets)!=len(refs):raise ValueError('invalid reference offsets')
    candidate_offset=_matrix([p.get('candidate_offset',[0]*n)],'candidate_offset',columns=n,high=q)[0]
    bm=p.get('block_map');target=integer(p.get('minimum_distance'),'minimum_distance')
    if not isinstance(bm,list) or len(bm)!=q or any(type(x) is not int or not 0<=x<1<<width for x in bm):raise ValueError('invalid map')
    block=np.asarray(bm,dtype=np.int64)
    witnesses=result.get('rejection_witnesses');candidates=result.get('candidates')
    if not isinstance(witnesses,list) or not isinstance(candidates,list):raise ValueError('invalid evidence lists')
    context.limit('max_trials',len(witnesses)+len(candidates));checked=0
    def vec(v,length,maximum):
        if not isinstance(v,list) or len(v)!=length or any(type(x) is not int or not 0<=x<maximum for x in v):raise ValueError('invalid evidence vector')
        return np.asarray(v,dtype=np.int64)
    for w in witnesses:
        if not isinstance(w,dict) or type(w.get('reference')) is not int or not 0<=w['reference']<len(refs):raise ValueError('invalid witness reference')
        a=vec(w.get('message_a'),k,q);b=vec(w.get('message_b'),k,q);t=vec(w.get('rejected_translate'),n,1<<width)
        d=_weight(t^refs[w['reference']]^_difference(a,b,G,q,block,candidate_offset,ring_offsets[w['reference']]))
        if type(w.get('distance')) is not int or d!=w['distance'] or d>=target:raise ValueError('invalid rejection witness')
        checked+=1
    exact=[c for c in candidates if isinstance(c,dict) and c.get('exhaustive_cross_check') is True]
    total=(q**k)**2*len(refs)*len(exact);context.limit('max_pairs',total)
    if exact:
        messages=np.array(list(product(range(q),repeat=k)),dtype=np.int64);images_a=block[(messages@G+candidate_offset)%q]
        for c in exact:
            t=vec(c.get('translate'),n,1<<width)
            for ri,ref in enumerate(refs):
                images_b=block[(messages@G+ring_offsets[ri])%q]
                for a in images_a:
                    for b in images_b:
                        if _weight(t^ref^a^b)<target:raise ValueError('false exhaustive cross-distance claim')
    return Verification('unverified','replayed rejection witnesses and exhaustive cross checks; parent proof is separate',
                        {'facts':{},'witnesses_checked':checked,'cross_checks_replayed':len(exact)})

OPERATORS={'search.binary_translates':binary_translate_search}
INPUT_KINDS={'search.binary_translates':()}
OUTPUT_KINDS={'search.binary_translates':('search_result',)}
