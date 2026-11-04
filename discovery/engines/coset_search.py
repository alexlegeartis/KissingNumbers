"""Actual candidate search for Gray images of parity-kernel ring coset unions.

This API deliberately returns search data, never a certified spherical artifact.
All sampled distances are upper estimates. A separate exact proof is required.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from collections import deque
from itertools import product, combinations
from pathlib import Path
import time
import numpy as np
from discovery.core import Artifact, Verification, RunContext, digest, canonical_json


@dataclass(frozen=True)
class CosetAssemblyRequest:
    generator: tuple[tuple[int, ...], ...]
    parity_checks: tuple[tuple[int, ...], ...]
    modulus: int
    target_distance: int
    required_cosets: int
    seed: int = 0
    offsets_per_family: int = 128
    information_sets: int = 8
    refinement_sets: int = 128
    max_vertices: int = 128
    max_cliques: int = 64
    max_corrections: int = 64
    companion_cardinality_forecast: int = 0
    families: tuple[str, ...] = ('uniform', 'even', 'quadratic')


@dataclass
class CosetAssemblySearch:
    status: str
    request: dict
    forecast: dict
    candidates: list[dict]
    rejection_witnesses: list[dict]
    statistics: dict
    torsion_guard: dict

    def to_dict(self):
        return asdict(self)


def _integer(x, name, minimum=1):
    if type(x) is not int or x < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return x


def _binary_rank(A):
    A=np.asarray(A,dtype=np.int64).copy()%2;at=0
    for j in range(A.shape[1]):
        ids=np.flatnonzero(A[at:,j])
        if not len(ids):continue
        i=at+int(ids[0]);A[[at,i]]=A[[i,at]]
        for k in range(len(A)):
            if k!=at and A[k,j]:A[k]^=A[at]
        at+=1
        if at==len(A):break
    return at


def _information_set(G,q,permutation):
    k,n=G.shape;A=G.copy();R=np.eye(k,dtype=np.int64);pivots=[];at=0
    for col in permutation:
        ids=np.flatnonzero(A[at:,col]%2)
        if not len(ids):continue
        j=at+int(ids[0]);A[[at,j]]=A[[j,at]];R[[at,j]]=R[[j,at]]
        inv=pow(int(A[at,col]),-1,q);A[at]=A[at]*inv%q;R[at]=R[at]*inv%q
        for i in range(k):
            if i!=at:
                c=int(A[i,col]);A[i]=(A[i]-c*A[at])%q;R[i]=(R[i]-c*R[at])%q
        pivots.append(int(col));at+=1
        if at==k:break
    if at!=k:raise ValueError('generator must have full residue row rank')
    return np.array(pivots),A,R


def _canonical(V,base,q):
    I,A,R=base
    return (V-V[:,I]@A)%q


def _dedup(V,L,base,q):
    I,A,R=base;keys=np.column_stack([_canonical(V,base,q),(V[:,I]@R@L.T)%2])
    _,indices=np.unique(keys,axis=0,return_index=True)
    return V[np.sort(indices)]


def _plans(G,L,q,count,rng,correction_limit):
    r=len(L);plans=[]
    for _ in range(count):
        I,A,R=_information_set(G,q,rng.permutation(G.shape[1]));U=L@R.T%2
        selected=[]
        for j in rng.permutation(G.shape[0]):
            trial=selected+[int(j)]
            if _binary_rank(U[:,trial].T)>len(selected):selected=trial
            if len(selected)==r:break
        codes=[sum(int(U[i,j])<<i for i in range(r))for j in selected]
        solve=np.zeros((1<<r,r),dtype=np.int64)
        for label in range(1<<r):
            syndrome=0
            for j,c in enumerate(codes):
                if label>>j&1:syndrome^=c
            solve[syndrome]=[(label>>j)&1 for j in range(r)]
        total=(q//2)**r
        if total<=correction_limit:labels=range(total)
        else:
            selected_labels={0}
            while len(selected_labels)<correction_limit:selected_labels.add(int(rng.integers(1,total)))
            labels=sorted(selected_labels)
        corrections=[[(int(label)//(q//2)**j)%(q//2)*2 for j in range(r)]for label in labels]
        plans.append((I,A,R,selected,solve,corrections))
    return plans


def _screen(V,L,q,plans,deadline=None):
    k=L.shape[1];r=len(L);best=np.full(len(V),V.shape[1]*q//2+1,dtype=np.int64);messages=np.zeros((len(V),k),dtype=np.int64)
    weight=np.full(q,q//4,dtype=np.int64);weight[0]=0;weight[q//2]=q//2
    for I,A,R,selected,solve,corrections in plans:
        if deadline is not None and time.monotonic()>deadline:break
        chosen=V[:,I];coeff=chosen@R%q;syndrome=coeff@L.T%2
        labels=sum(syndrome[:,j]*(1<<j)for j in range(r));bits=solve[labels];residual=(V-chosen@A)%q
        for correction in corrections:
            delta=bits+np.array(correction);W=(residual-delta@A[selected])%q;distance=weight[W].sum(axis=1);improve=distance<best
            if np.any(improve):best[improve]=distance[improve];messages[improve]=(coeff[improve]+delta[improve]@R[selected])%q
    if np.any(messages@L.T%2):raise ValueError('constrained oracle membership failure')
    return best,messages


class _SearchCutoff(RuntimeError):
    pass


def _clique(adjacency,target,deadline=None,max_branches=None):
    branches=0
    def visit(P,chosen):
        nonlocal branches
        branches+=1
        if (deadline is not None and time.monotonic()>deadline)or(max_branches is not None and branches>max_branches):raise _SearchCutoff('clique branch budget exhausted')
        if len(chosen)==target:return chosen
        while P:
            if P.bit_count()<target-len(chosen):return None
            bit=P&-P;v=bit.bit_length()-1;P^=bit
            result=visit(P&adjacency[v],chosen+[v])
            if result is not None:return result
        return None
    return visit((1<<len(adjacency))-1,[])


def _covering_radius(G,base,q,limit,deadline):
    k,n=G.shape;outside=[j for j in range(n)if j not in base[0]]
    if len(outside)>16:return None
    states=1<<len(outside)
    if states>min(limit,65536):return None
    columns=[]
    for j in range(n):
        e=np.zeros((1,n),dtype=np.int64);e[0,j]=1;s=_canonical(e,base,q)[0,outside]%2
        columns.append(sum(int(v)<<i for i,v in enumerate(s)))
    D=[-1]*states;D[0]=0;Q=deque([0]);checked=0
    while Q:
        a=Q.popleft()
        for b in columns:
            if D[a^b]<0:D[a^b]=D[a]+1;Q.append(a^b)
        checked+=1
        if checked%1024==0 and time.monotonic()>deadline:return None
    return max(D)


def _validate_request(request,context):
    if not isinstance(request,CosetAssemblyRequest):raise ValueError('typed CosetAssemblyRequest required')
    q=request.modulus
    if type(q)is not int or q not in(4,8,16):raise ValueError('modulus must be4,8,or16')
    for name in ('target_distance','required_cosets','offsets_per_family','information_sets','refinement_sets','max_vertices','max_cliques','max_corrections'):_integer(getattr(request,name),name)
    _integer(request.seed,'seed',0);_integer(request.companion_cardinality_forecast,'companion_cardinality_forecast',0)
    raw=request.generator
    if not isinstance(raw,(tuple,list))or not raw or not isinstance(raw[0],(tuple,list))or not raw[0]:raise ValueError('nonempty generator matrix required')
    k,n=len(raw),len(raw[0]);context.limit('max_dimension',n*q//2)
    if any(not isinstance(row,(tuple,list))or len(row)!=n or any(type(x)is not int or not 0<=x<q for x in row)for row in raw):raise ValueError('generator dimensions/domain')
    checks=request.parity_checks
    if not isinstance(checks,(tuple,list)):raise ValueError('parity matrix required')
    r=len(checks)
    if not 1<=r<=min(k,8)or any(not isinstance(row,(tuple,list))or len(row)!=k or any(type(x)is not int or x not in(0,1)for x in row)for row in checks):raise ValueError('binary parity matrix dimensions/domain')
    context.limit('max_rational_bits',k*(q.bit_length()-1)+2);context.limit('max_trials',request.max_cliques)
    context.limit('max_trials',request.max_corrections)
    context.limit('max_trials',request.information_sets+request.refinement_sets)
    context.limit('max_gram_entries',(k*n+k*k)*(request.information_sets+request.refinement_sets))
    if (1<<r)-1>request.max_vertices:raise ValueError('graph budget must include all parity-quotient controls')
    context.limit('max_pairs',request.max_vertices*(request.max_vertices-1)//2)
    context.limit('max_gram_entries',(k*n+k*k)*(request.information_sets+request.refinement_sets)+request.max_vertices*(request.max_vertices-1)//2*n)
    if request.target_distance>n*q//2:raise ValueError('target distance exceeds Gray length')
    if not isinstance(request.families,(tuple,list))or any(not isinstance(x,str)for x in request.families):raise ValueError('offset family names required')
    context.limit('max_rational_bits',request.companion_cardinality_forecast.bit_length())
    context.limit('max_codewords',len(request.families)*(request.offsets_per_family+1)*(1<<r))
    if request.required_cosets<2 or request.required_cosets>request.max_vertices+1:raise ValueError('required_cosets outside graph budget')
    if not request.families or len(set(request.families))!=len(request.families)or set(request.families)-{'uniform','even','quadratic','torsion'}:raise ValueError('unknown or duplicate offset family')
    G=np.array(raw,dtype=np.int64);L=np.array(checks,dtype=np.int64)
    if _binary_rank(L)!=r:raise ValueError('parity checks must be independent')
    base=_information_set(G,q,range(n))
    return G,L,q,k,n,r,base


def search_gray_coset_assemblies(request: CosetAssemblyRequest, context=None) -> CosetAssemblySearch:
    """Generate and refine real compressed candidate assemblies.

    A candidate consists of cosets of D={message*G: L*message=0 mod2}.
    It carries no distance lower bound, no spherical-code facts, and no novelty
    claim. Companion cardinality is an explicit forecast assumption only.
    """
    context=context or RunContext(Path('.'));start=time.monotonic();deadline=start+context.budget['solver_seconds']
    G,L,q,k,n,r,base=_validate_request(request,context)
    outside=[j for j in range(n)if j not in base[0]];rng=np.random.default_rng(request.seed)
    radius=_covering_radius(G,base,q,context.budget['max_codewords'],deadline)
    guard_active=radius is not None and (q//2)*radius<request.target_distance and request.required_cosets>(1<<r)
    guard={'binary_covering_radius':radius,'active':guard_active,'scope':'pure half-modulus-torsion syndromes only','maximum_compatible_cosets_in_that_family':1<<r,'proof':'same parity-label differences have homogeneous distance<=modulus/2*binary_covering_radius; color by kernel quotient label'}
    prototypes=[np.zeros((1,n),dtype=np.int64)]
    bits=((np.arange(n)[:,None]>>np.arange(max(1,(n-1).bit_length())))&1).T
    F=np.vstack([np.ones(n,dtype=int),bits,*[bits[i]*bits[j]for i in range(len(bits))for j in range(i)]])
    for family in request.families:
        V=np.zeros((request.offsets_per_family,n),dtype=np.int64)
        if family=='quadratic':V=_canonical(rng.integers(0,q,(len(V),len(F)))@F%q,base,q)
        else:
            scale=2 if family=='even'else q//2 if family=='torsion'else 1
            V[:,outside]=scale*rng.integers(0,q//scale,(len(V),len(outside)))
        prototypes.append(V)
    prototypes=np.unique(np.vstack(prototypes),axis=0);removed=0
    if guard_active and set(request.families)=={'torsion'}:
        discard=np.all(prototypes%(q//2)==0,axis=1)&np.any(prototypes,axis=1);removed=int(discard.sum());prototypes=prototypes[~discard]
    guard['prototypes_removed']=removed
    guard['mixed_requests_retain_torsion_vertices']=True
    # Message-coordinate columns spanning the parity quotient provide all labels.
    selected=[]
    for j in range(k):
        if _binary_rank(L[:,selected+[j]].T)>len(selected):selected.append(j)
        if len(selected)==r:break
    lifts=[sum((G[selected[j]]for j in range(r)if label>>j&1),np.zeros(n,dtype=np.int64))%q for label in range(1<<r)]
    V=_dedup(np.vstack([(prototypes+e)%q for e in lifts]),L,base,q);V=V[np.any(V,axis=1)]
    plans=_plans(G,L,q,request.information_sets+request.refinement_sets,rng,request.max_corrections)
    d,msg=_screen(V,L,q,plans[:request.information_sets],deadline);eligible=np.flatnonzero(d>=request.target_distance)
    controls=[int(np.flatnonzero(np.all(V==e,axis=1))[0])for e in lifts[1:]]
    chosen=sorted(set(controls+list(eligible[np.argsort(d[eligible])[-request.max_vertices:]])))
    chosen=[i for i in chosen if d[i]>=request.target_distance]
    if len(chosen)>request.max_vertices:
        rest=[i for i in chosen if i not in controls];chosen=controls+rest[:request.max_vertices-len(controls)]
    S=V[chosen];pairs=list(combinations(range(len(S)),2));context.limit('max_pairs',len(pairs))
    DV=np.array([(S[i]-S[j])%q for i,j in pairs],dtype=np.int64)
    pd,pm=_screen(DV,L,q,plans[:request.information_sets],deadline)if pairs else(np.array([],dtype=int),np.zeros((0,k),dtype=int))
    adj=[0]*len(S);torsion_edges_pruned=0
    pure_torsion=np.all(_canonical(S,base,q)%(q//2)==0,axis=1)
    labels=(S[:,base[0]]@base[2]@L.T)%2
    for (i,j),distance in zip(pairs,pd):
        same_torsion_color=guard_active and pure_torsion[i]and pure_torsion[j]and np.array_equal(labels[i],labels[j])
        if same_torsion_color:torsion_edges_pruned+=1;continue
        if distance>=request.target_distance:adj[i]|=1<<j;adj[j]|=1<<i
    guard['same_color_torsion_edges_pruned']=torsion_edges_pruned
    candidates=[];witnesses=[];attempts=0;status='unverified_search_complete'
    for attempt in range(request.max_cliques):
        if time.monotonic()>deadline:status='unverified_budget_cutoff';break
        try:clique=_clique(adj,request.required_cosets-1,deadline,context.budget['max_rank_work'])
        except _SearchCutoff:status='unverified_budget_cutoff';break
        if clique is None:break
        attempts+=1;offsets=np.vstack([np.zeros((1,n),dtype=np.int64),S[clique]]);edges=list(combinations(range(len(offsets)),2));diff=np.array([(offsets[i]-offsets[j])%q for i,j in edges])
        strong,sm=_screen(diff,L,q,plans[request.information_sets:],deadline);bad=np.flatnonzero(strong<request.target_distance)
        if not len(bad):
            candidates.append({'offsets':offsets.tolist(),'parity_checks':L.tolist(),'distance_status':'UNVERIFIED','sampled_minimum_upper_estimates':strong.tolist(),'exact_proof_obligations':['internal parent/subcode distance','all inter-coset minimum distances','coset distinctness','companion compatibility if used']});break
        for index in bad:
            i,j=edges[int(index)];distance=int(strong[index])
            if distance<=0:raise ValueError('quotient deduplication invariant failed')
            witnesses.append({'offset_difference':diff[index].tolist(),'parent_message':sm[index].tolist(),'distance':distance})
            if i==0:
                v=clique[j-1]
                for a in range(len(S)):adj[a]&=~(1<<v)
                adj[v]=0
            else:
                a,b=clique[i-1],clique[j-1];adj[a]&=~(1<<b);adj[b]&=~(1<<a)
    else:status='unverified_clique_budget_cutoff'
    parent=q**k;size=parent//(1<<r)
    forecast={'parent_cardinality':parent,'subcode_index':1<<r,'subcode_cardinality':size,'required_cosets':request.required_cosets,'union_cardinality_if_distinct':size*request.required_cosets,'companion_cardinality_assumption':request.companion_cardinality_forecast,'total_cardinality_if_certified':size*request.required_cosets+request.companion_cardinality_forecast,'gray_binary_length':n*q//2,'scope':'cardinality forecast only; not verified spherical facts'}
    statistics={'quotient_distinct_offsets':len(V),'eligible_upper_estimates':len(eligible),'retained_vertices':len(S),'pair_differences_queried':len(pairs),'cliques_refined':attempts,'elapsed_seconds':time.monotonic()-start,'request_sha256':digest(asdict(request))}
    return CosetAssemblySearch(status,asdict(request),forecast,candidates,witnesses,statistics,guard)


BACKEND='gray_coset_assembly'
FORMAT='gray_coset_assembly_search_v1'
_PROOF_OBLIGATIONS=['internal parent/subcode distance','all inter-coset minimum distances',
                    'coset distinctness','companion compatibility if used']


def request_from_params(params):
    if not isinstance(params,dict)or set(params)-set(CosetAssemblyRequest.__dataclass_fields__):
        raise ValueError('unknown coset assembly request fields')
    converted=dict(params)
    try:
        for name in ('generator','parity_checks'):
            if name in converted:converted[name]=tuple(tuple(row)for row in converted[name])
        if 'families'in converted:converted['families']=tuple(converted['families'])
        return CosetAssemblyRequest(**converted)
    except (TypeError,ValueError)as exc:
        raise ValueError('invalid compressed coset request')from exc


def coset_assembly(inputs,params,context):
    if inputs:raise ValueError('coset_assembly takes compressed parameters, no input artifacts')
    result=search_gray_coset_assemblies(request_from_params(params),context)
    return Artifact('search_result',{'backend':BACKEND,'format':FORMAT,'result':result.to_dict()},
        {},Verification('unverified','bounded coset candidate generation; exact distance proof absent',
                        {'facts':{},'cardinality_forecast_is_not_certified':True}))


def verify_artifact(artifact,context=None):
    """Validate search schema and negative witnesses, always return UNVERIFIED.

    Stored forecasts, candidate estimates and status never become sphere facts.
    This is deliberately not an exact positive-distance verifier.
    """
    context=context or RunContext(Path('.'));data=artifact.data
    if artifact.kind!='search_result'or artifact.facts!={}:
        raise ValueError('coset search has no certified geometric facts')
    if set(data)!={'backend','format','result'}or data['backend']!=BACKEND or data['format']!=FORMAT:
        raise ValueError('invalid coset search artifact schema')
    result=data['result'];expected=set(CosetAssemblySearch.__dataclass_fields__)
    if not isinstance(result,dict)or set(result)!=expected:raise ValueError('invalid search result fields')
    if result['status']not in('unverified_search_complete','unverified_budget_cutoff','unverified_clique_budget_cutoff'):
        raise ValueError('coset search status cannot assert certification')
    request=request_from_params(result['request']);G,L,q,k,n,r,base=_validate_request(request,context)
    if not isinstance(result['forecast'],dict):raise ValueError('forecast object required')
    if result['forecast'].get('scope')!='cardinality forecast only; not verified spherical facts':
        raise ValueError('forecast must remain explicitly conditional')
    forecast_fields={'parent_cardinality','subcode_index','subcode_cardinality','required_cosets',
        'union_cardinality_if_distinct','companion_cardinality_assumption','total_cardinality_if_certified',
        'gray_binary_length','scope'}
    if set(result['forecast'])!=forecast_fields:raise ValueError('unknown forecast fields')
    for key,value in result['forecast'].items():
        if key!='scope':_integer(value,'forecast '+key,0)
    parent=q**k;size=parent//(1<<r)
    expected_values={'parent_cardinality':parent,'subcode_index':1<<r,'subcode_cardinality':size,
        'required_cosets':request.required_cosets,'union_cardinality_if_distinct':size*request.required_cosets,
        'companion_cardinality_assumption':request.companion_cardinality_forecast,
        'total_cardinality_if_certified':size*request.required_cosets+request.companion_cardinality_forecast,
        'gray_binary_length':n*q//2}
    if any(result['forecast'][key]!=value for key,value in expected_values.items()):
        raise ValueError('conditional forecast inconsistent with request')
    def matrix(rows,width,bound,name):
        if not isinstance(rows,list)or any(not isinstance(row,list)or len(row)!=width or
            any(type(x)is not int or not 0<=x<bound for x in row)for row in rows):
            raise ValueError(name+' dimensions/domain')
    candidates=result['candidates']
    if not isinstance(candidates,list)or len(candidates)>1:raise ValueError('invalid bounded candidate list')
    for candidate in candidates:
        if not isinstance(candidate,dict)or set(candidate)!={'offsets','parity_checks','distance_status','sampled_minimum_upper_estimates','exact_proof_obligations'}:
            raise ValueError('unknown candidate fields')
        if candidate['distance_status']!='UNVERIFIED'or candidate['exact_proof_obligations']!=_PROOF_OBLIGATIONS:
            raise ValueError('candidate proof obligations cannot be removed')
        matrix(candidate['offsets'],n,q,'candidate offsets')
        if len(candidate['offsets'])!=request.required_cosets or candidate['parity_checks']!=L.tolist():
            raise ValueError('candidate construction does not match request')
        offsets=np.array(candidate['offsets'],dtype=np.int64)
        if len(_dedup(offsets,L,base,q))!=len(offsets):raise ValueError('candidate has coincident subcode cosets')
        estimates=candidate['sampled_minimum_upper_estimates']
        if not isinstance(estimates,list)or len(estimates)!=len(offsets)*(len(offsets)-1)//2 or any(
            type(x)is not int or not request.target_distance<=x<=n*q//2+1 for x in estimates):
            raise ValueError('invalid sampled upper estimates')
    witnesses=result['rejection_witnesses']
    if not isinstance(witnesses,list):raise ValueError('rejection witness list required')
    context.limit('max_pairs',len(witnesses))
    for witness in witnesses:
        if not isinstance(witness,dict)or set(witness)!={'offset_difference','parent_message','distance'}:raise ValueError('unknown rejection fields')
        matrix([witness['offset_difference']],n,q,'witness offset');matrix([witness['parent_message']],k,q,'witness message')
        c=np.array(witness['parent_message'],dtype=np.int64)
        if np.any(L@c%2):raise ValueError('rejection message outside subcode')
        word=(np.array(witness['offset_difference'])-c@G)%q
        actual=sum(0 if int(x)==0 else q//2 if int(x)==q//2 else q//4 for x in word)
        if type(witness['distance'])is not int or not 0<actual==witness['distance']<request.target_distance:
            raise ValueError('invalid actual rejection witness')
    statistics_fields={'quotient_distinct_offsets','eligible_upper_estimates','retained_vertices','pair_differences_queried','cliques_refined','elapsed_seconds','request_sha256'}
    if not isinstance(result['statistics'],dict)or set(result['statistics'])!=statistics_fields or result['statistics'].get('request_sha256')!=digest(asdict(request)):
        raise ValueError('search statistics are not bound to request')
    for key,value in result['statistics'].items():
        if key not in('elapsed_seconds','request_sha256'):_integer(value,'statistics '+key,0)
    elapsed=result['statistics']['elapsed_seconds']
    if type(elapsed)not in(int,float)or not 0<=elapsed<float('inf'):raise ValueError('invalid elapsed time')
    guard_fields={'binary_covering_radius','active','scope','maximum_compatible_cosets_in_that_family','proof','prototypes_removed','mixed_requests_retain_torsion_vertices','same_color_torsion_edges_pruned'}
    if not isinstance(result['torsion_guard'],dict)or set(result['torsion_guard'])!=guard_fields:raise ValueError('torsion scope metadata required')
    guard=result['torsion_guard']
    if guard['scope']!='pure half-modulus-torsion syndromes only'or type(guard['active'])is not bool or guard['mixed_requests_retain_torsion_vertices']is not True:raise ValueError('invalid torsion scope')
    if guard['binary_covering_radius']is not None:_integer(guard['binary_covering_radius'],'covering radius',0)
    if guard['maximum_compatible_cosets_in_that_family']!=1<<r:raise ValueError('wrong torsion family index')
    for key in('prototypes_removed','same_color_torsion_edges_pruned'):_integer(guard[key],key,0)
    if set(request.families)!={'torsion'}and result['torsion_guard'].get('prototypes_removed')!=0:
        raise ValueError('mixed assemblies must retain torsion vertices')
    return Verification('unverified','search schema and rejection witnesses checked; positive distance proof absent',
        {'facts':{},'candidate_count':len(candidates),'cardinality_forecast_is_not_certified':True,
         'proof_obligations_remain':_PROOF_OBLIGATIONS})


OPERATORS={'search.coset_assembly':coset_assembly}
INPUT_KINDS={'search.coset_assembly':()}
OUTPUT_KINDS={'search.coset_assembly':('search_result',)}
