"""Synthesize orthogonal representations from finite-group coset actions.

Search the new representation using the ordinary joint group-orbit optimizer.
The output remains provisional until the separate exact orbit promotion.
"""
import math
from collections import deque
import numpy as np
from ..core import Artifact,RunContext,Verification,check_params,integer
from .group_seed_search import _closure,search_group_seeds
from .exact_orbits import _generators


def coset_generators(generators,subgroup_generators,dimension,copies,trivial,context,supplemental_generators=None,supplemental_dimension=0):
    G,den=_closure(generators,dimension,context)
    H,subden=_closure(subgroup_generators,dimension,context)
    if subden!=den:raise ValueError('group and subgroup denominators must agree')
    keys={tuple(m.ravel()):i for i,m in enumerate(G)}
    if any(tuple(h.ravel()) not in keys for h in H):raise ValueError('subgroup is not contained in group')
    owner={};reps=[]
    for i,g in enumerate(G):
        if i in owner:continue
        q=len(reps); reps.append(i)
        for h in H:
            p=g@h
            if np.any(p%den):raise ValueError('group lattice is not closed')
            k=keys.get(tuple((p//den).ravel()))
            if k is None or k in owner:raise ValueError('invalid coset partition')
            owner[k]=q
    m=len(reps);context.limit('max_dimension',m)
    if len(owner)!=len(G) or m*len(H)!=len(G):raise ValueError('invalid subgroup index')
    k=math.isqrt(m)
    if k*k!=m or k<2:raise ValueError('standard rational representation requires a nontrivial square index')
    hden=k*(k-1);w=np.ones(m,dtype=np.int64);w[-1]=1-k;T=hden*np.eye(m,dtype=np.int64)-np.outer(w,w)
    # A second representation can be glued generator-by-generator only after
    # checking all relations of the source group on it, via synchronized BFS.
    extra=[];extra_den=1
    if supplemental_generators is not None:
        ed=integer(supplemental_dimension,'supplemental_dimension');context.limit('max_dimension',ed)
        ex=_generators(supplemental_generators,ed,context)
        if len(ex)!=len(generators):raise ValueError('supplemental generators must align with source generators')
        extra_den=ex[0][1]
        if any(q!=extra_den for _,q in ex):raise ValueError('supplemental generators must share a denominator')
        if ed*extra_den*extra_den>=2**62:raise ValueError('supplemental products exceed safe range')
        extra=[np.asarray(x,dtype=np.int64) for x,_ in ex]
        source=[np.asarray(x,dtype=np.int64) for x,_ in _generators(generators,dimension,context)]
        identity=np.eye(dimension,dtype=np.int64)*den;ei=np.eye(ed,dtype=np.int64)*extra_den
        seen={tuple(identity.ravel()):tuple(ei.ravel())};queue=deque([(identity,ei)])
        while queue:
            a,b=queue.popleft()
            for g,e in zip(source,extra):
                aa=g@a;bb=e@b
                if np.any(aa%den) or np.any(bb%extra_den):raise ValueError('supplemental representation leaves its rational lattice')
                aa//=den;bb//=extra_den;ak=tuple(aa.ravel());bk=tuple(bb.ravel())
                if ak in seen:
                    if seen[ak]!=bk:raise ValueError('supplemental generators violate source-group relations')
                else:seen[ak]=bk;queue.append((aa,bb))
        if len(seen)!=len(G):raise ValueError('incomplete representation relation check')
    elif supplemental_dimension!=0:raise ValueError('supplemental dimension requires generators')
    outdim=(m-1)*copies+trivial+supplemental_dimension;context.limit('max_dimension',outdim)
    if m*hden*hden>=2**60:raise ValueError('coset basis arithmetic exceeds safe range')
    raw=_generators(generators,dimension,context);result=[];qden=hden*hden;common=math.lcm(qden,extra_den)
    if outdim*common*common>=2**62:raise ValueError('combined representation exceeds safe range')
    for n,(mat,d) in enumerate(raw):
        g=np.asarray(mat,dtype=np.int64); perm=[]
        for i in reps:
            prod=g@G[i]
            if np.any(prod%den):raise ValueError('group lattice is not closed')
            perm.append(owner[keys[tuple((prod//den).ravel())]])
        if len(set(perm))!=m:raise ValueError('coset action is not a permutation')
        P=np.eye(m,dtype=np.int64)[:,perm];Q=T@P@T
        if np.any(Q[:-1,-1]) or np.any(Q[-1,:-1]) or Q[-1,-1]!=qden:raise ValueError('standard subspace is not preserved')
        R=Q[:-1,:-1]; M=np.zeros((outdim,outdim),dtype=np.int64)
        for c in range(copies):M[c*(m-1):(c+1)*(m-1),c*(m-1):(c+1)*(m-1)]=R*(common//qden)
        begin=(m-1)*copies
        if extra:M[begin:begin+supplemental_dimension,begin:begin+supplemental_dimension]=extra[n]*(common//extra_den)
        for j in range(begin+supplemental_dimension,outdim):M[j,j]=common
        result.append({'matrix':M.tolist(),'denominator':common})
    return result,{'group_order':len(G),'subgroup_order':len(H),'index':m,'standard_dimension':m-1,'representation_dimension':outdim}


def search_coset_orbits(inputs,params,context):
    if inputs:raise ValueError('coset search takes no inputs')
    allowed=('generators','subgroup_generators','dimension','copies','trivial_dimensions','orbit_count','attempts','max_iterations','seed','seconds','scale','initial_seeds','supplemental_generators','supplemental_dimension')
    check_params(params,allowed)
    d=integer(params.get('dimension'),'dimension');copies=integer(params.get('copies',1),'copies');trivial=integer(params.get('trivial_dimensions',0),'trivial_dimensions',0)
    extra_d=integer(params.get('supplemental_dimension',0),'supplemental_dimension',0)
    gens,info=coset_generators(params.get('generators'),params.get('subgroup_generators'),d,copies,trivial,context,params.get('supplemental_generators'),extra_d)
    p={k:params[k] for k in ('orbit_count','attempts','max_iterations','seed','seconds','scale','initial_seeds') if k in params};p['generators']=gens;p['dimension']=info['representation_dimension']
    result=search_group_seeds(p,context);result['coset_representation']=info
    return Artifact('search_result',{'backend':'group_seed_search','format':'group_seed_search_v1','result':result},{},Verification('unverified','coset-representation search proposal; exact promotion required',{}))

OPERATORS={'search.coset_group_orbits':search_coset_orbits}
INPUT_KINDS={'search.coset_group_orbits':()}
OUTPUT_KINDS={'search.coset_group_orbits':('search_result',)}
