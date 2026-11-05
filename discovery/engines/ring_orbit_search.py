"""Search partial 2C-orbit assemblies for affine Z8 code isometries.

Witnessed conflicts generate a bipartite graph between C/2C and (αC+t)/2C.
A matching is an exact upper bound on any whole-orbit assembly, even when the
oracle has missed conflicts. Independent sets in the sampled graph remain
provisional until every omitted cross-distance has been proved safe.
"""
from collections import deque
from pathlib import Path
import math,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from ..core import Artifact,Verification,RunContext,check_params,integer
from .binary_translate_search import _matrix


def _parameters(p,context):
    check_params(p,{'generator','alpha','shift','minimum_distance','max_orbits','attempts','steps','seed','seconds','initial_witnesses'})
    G=_matrix(p.get('generator'),'generator',high=8);k,n=G.shape
    if k>16 or k>n:raise ValueError('generator rank exceeds supported range')
    context.limit('max_dimension',4*n);N=1<<k;context.limit('max_codewords',N)
    alpha=p.get('alpha');shift=p.get('shift')
    if not isinstance(alpha,list) or len(alpha)!=n or any(type(x) is not int or x not in (1,5) for x in alpha):raise ValueError('alpha entries must be 1 or 5')
    if not isinstance(shift,list) or len(shift)!=n or any(type(x) is not int or not 0<=x<8 for x in shift):raise ValueError('shift entries must be in Z8')
    target=integer(p.get('minimum_distance'),'minimum_distance')
    if target>4*n:raise ValueError('distance exceeds binary length')
    bits=(np.arange(N,dtype=np.int64)[:,None]>>np.arange(k))&1;binary=(bits@G)%2
    masks=[sum(int(x)<<j for j,x in enumerate(row)) for row in binary]
    if len(set(masks))!=N:raise ValueError('generator must have independent rows modulo 2')
    lookup={x:i for i,x in enumerate(masks)};A=sum((x==5)<<j for j,x in enumerate(alpha));K=[(h,lookup[A&w]) for h,w in enumerate(masks) if A&w in lookup]
    return G,np.asarray(alpha,np.int64),np.asarray(shift,np.int64),target,bits,K


def _expand(witnesses,G,alpha,shift,target,bits,K,context):
    k,n=G.shape;N=len(bits);mask=sum((int(x)==5)<<j for j,x in enumerate(alpha));residue=(bits@G)%2
    lookup={sum(int(x)<<j for j,x in enumerate(row)):i for i,row in enumerate(residue)}
    context.limit('max_pairs',len(witnesses)*len(K));wt=np.array([0,2,2,2,4,2,2,2],np.int64);edges=set()
    for item in witnesses:
        if not isinstance(item,dict) or set(item)!={'candidate','source_coefficients','distance'}:raise ValueError('invalid orbit witness')
        u=item['candidate'];b=item['source_coefficients'];distance=item['distance']
        if type(u) is not int or not 0<=u<N or not isinstance(b,list) or len(b)!=k or any(type(x) is not int or not 0<=x<8 for x in b) or type(distance) is not int:raise ValueError('invalid orbit witness values')
        b=np.asarray(b,np.int64);src=sum(int(b[i]&1)<<i for i in range(k));old=bits[u]
        for h,z in K:
            assert lookup[mask & sum(int(x)<<j for j,x in enumerate(residue[h]))]==z
            hb=bits[h];zb=bits[z];u2=u^h;b2=(b+hb+4*zb-2*(old*hb))%8;v2=sum(int(b2[i]&1)<<i for i in range(k))
            x=(alpha*(bits[u2]@G)+shift)%8;y=b2@G%8;d=int(wt[(x-y)%8].sum())
            if d!=distance or d>=target or v2!=(src^h):raise ValueError('invalid expanded orbit conflict witness')
            edges.add((u2,v2))
    return edges


def _independent_set(M,match):
    N=M.shape[0];reverse=np.full(N,-1,np.int64)
    for i,j in enumerate(match):
        if j>=0:reverse[j]=i
    rows=np.zeros(N,bool);cols=np.zeros(N,bool);queue=deque(np.flatnonzero(match<0).tolist());rows[list(queue)]=True
    while queue:
        i=queue.popleft()
        for j in M.indices[M.indptr[i]:M.indptr[i+1]]:
            if match[i]==j or cols[j]:continue
            cols[j]=True;r=reverse[j]
            if r>=0 and not rows[r]:rows[r]=True;queue.append(int(r))
    return np.flatnonzero(rows).tolist(),np.flatnonzero(~cols).tolist()


def search_ring_orbits(p,context=None):
    context=context or RunContext(Path('.'));G,alpha,shift,target,bits,K=_parameters(p,context);k,n=G.shape;N=len(bits)
    max_orbits=integer(p.get('max_orbits',min(64,N)),'max_orbits');context.limit('max_codewords',max_orbits)
    if max_orbits>N:raise ValueError('max_orbits exceeds total orbits')
    attempts=integer(p.get('attempts',2),'attempts');context.limit('max_trials',attempts)
    steps=integer(p.get('steps',500),'steps');context.limit('max_rank_work',max_orbits*attempts*steps*n*8)
    seed=integer(p.get('seed',0),'seed',0);seconds=p.get('seconds',10)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);rng=np.random.default_rng(seed);deadline=time.monotonic()+seconds;wt=np.array([0,2,2,2,4,2,2,2],np.int64)
    initial=p.get('initial_witnesses',[])
    if not isinstance(initial,list):raise ValueError('initial_witnesses must be a list')
    context.limit('max_codewords',len(initial));_expand(initial,G,alpha,shift,target,bits,K,context)
    order=rng.choice(N,max_orbits,replace=False);witnesses=[dict(x) for x in initial];sampled=0
    triangular=k<=n and not np.any(np.tril(G[:,:k],-1)) and all(math.gcd(int(G[i,i]),8)==1 for i in range(k))
    for u0 in order:
        if time.monotonic()>=deadline:break
        u=int(u0);x=(alpha*(bits[u]@G)+shift)%8;sampled+=1;found=None
        # Exact information-set decoding finds coincident points cheaply.
        if triangular:
            b=np.zeros(k,np.int64)
            for j in range(k):b[j]=((int(x[j])-int(b@G[:,j]))*pow(int(G[j,j]),-1,8))%8
            if np.array_equal(b@G%8,x):found=(b,0)
        for attempt in range(attempts):
            if found is not None or time.monotonic()>=deadline:break
            b=rng.integers(8,size=k);v=b@G%8
            for it in range(steps):
                d=int(wt[(x-v)%8].sum())
                if d<target:found=(b.copy(),d);break
                i=int(rng.integers(k));old=int(b[i]);choices=(v[None,:]+(np.arange(8)[:,None]-old)*G[i])%8;scores=wt[(x[None,:]-choices)%8].sum(axis=1)
                temperature=.4+2.6*(1-(it%1000)/1000);prob=np.exp(-(scores-scores.min())/temperature);prob/=prob.sum();j=int(rng.choice(8,p=prob));b[i]=j;v=choices[j]
        if found is not None:
            b,d=found;witnesses.append({'candidate':u,'source_coefficients':b.tolist(),'distance':int(d)})
    edges=_expand(witnesses,G,alpha,shift,target,bits,K,context)
    if edges:
        pairs=np.asarray(sorted(edges),dtype=np.int64);M=coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(N,N)).tocsr()
    else:M=coo_matrix((N,N)).tocsr()
    match=maximum_bipartite_matching(M,perm_type='column');matching=[[int(i),int(j)] for i,j in enumerate(match) if j>=0]
    selected_candidate,selected_source=_independent_set(M,match)
    return {'parameters':p,'sampled_candidate_orbits':sampled,'subgroup_symmetries':len(K),'orbit_size':4**k,
            'witnesses':witnesses,'matching':matching,'whole_orbit_count_upper_bound':2*N-len(matching),
            'provisional_candidate_orbits':selected_candidate,'provisional_source_orbits':selected_source,
            'scope':'Matching witnesses prove an upper bound for this whole-orbit pool. The sampled-graph independent set is provisional; unobserved conflicts and parent-code validity require separate proof.'}


def ring_orbit_search(inputs,params,context):
    if inputs:raise ValueError('ring orbit search takes no inputs')
    result=search_ring_orbits(params,context)
    return Artifact('search_result',{'backend':'ring_orbit_search','format':'ring_orbit_search_v1','result':result},{},Verification('unverified','orbit assembly search; parent code and missing-edge obligations remain',{}))


def verify_artifact(artifact,context=None):
    context=context or RunContext(Path('.'))
    if artifact.kind!='search_result' or artifact.data.get('backend')!='ring_orbit_search' or artifact.data.get('format')!='ring_orbit_search_v1':raise ValueError('invalid orbit search receipt')
    r=artifact.data.get('result');p=r.get('parameters') if isinstance(r,dict) else None
    if not isinstance(p,dict):raise ValueError('missing parameters')
    G,alpha,shift,target,bits,K=_parameters(p,context);N=len(bits);W=r.get('witnesses');matching=r.get('matching')
    if not isinstance(W,list) or not isinstance(matching,list):raise ValueError('invalid orbit evidence')
    context.limit('max_codewords',len(W));edges=_expand(W,G,alpha,shift,target,bits,K,context);used_a=set();used_b=set()
    for pair in matching:
        if not isinstance(pair,list) or len(pair)!=2 or any(type(x) is not int or not 0<=x<N for x in pair):raise ValueError('invalid matching pair')
        a,b=pair
        if (a,b) not in edges or a in used_a or b in used_b:raise ValueError('matching is not supported by disjoint witnessed edges')
        used_a.add(a);used_b.add(b)
    if r.get('whole_orbit_count_upper_bound')!=2*N-len(matching) or r.get('orbit_size')!=4**G.shape[0]:raise ValueError('invalid orbit bound')
    return Verification('unverified','replayed orbit conflict witnesses and disjoint matching; missing-edge and parent-code proofs remain',
                        {'facts':{},'witnesses_checked':len(W),'matching_edges_checked':len(matching)})

OPERATORS={'search.ring_orbit_assemblies':ring_orbit_search}
INPUT_KINDS={'search.ring_orbit_assemblies':()}
OUTPUT_KINDS={'search.ring_orbit_assemblies':('search_result',)}
