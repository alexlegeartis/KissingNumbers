"""Finite-field polynomial images and partial signed-octad exchange search.

Neither a support pool nor an exchange proposal is a spherical-code proof.
The exchange backend checks its selected lines against the supplied finite
source description; certification of the complete source is separate.
"""
from itertools import combinations,product
from math import comb,isfinite
import numpy as np
from ..core import Artifact,Verification,check_params,integer,digest


def _mul(a,b,m,mod):
    q=1<<m;out=0
    while b:
        if b&1:out^=a
        a<<=1
        if a&q:a^=mod
        b>>=1
    return out


def _power(x,n,m,mod):
    out=1
    while n:
        if n&1:out=_mul(out,x,m,mod)
        x=_mul(x,x,m,mod);n>>=1
    return out


def _field(m,mod):
    if type(m) is not int or not 3<=m<=6:raise ValueError('m must be 3 through 6')
    q=1<<m
    if type(mod) is not int or not(q<=mod<2*q) or not(mod&1):raise ValueError('modulus must be an odd monic degree-m binary polynomial')
    # Exhaustive nonzero multiplication checks irreducibility for these small fields.
    if any(len({_mul(a,b,m,mod) for b in range(1,q)})!=q-1 or any(_mul(a,b,m,mod)==0 for b in range(1,q)) for a in range(1,q)):
        raise ValueError('reducible modulus')
    return q


def affine_flats(m,k,context):
    q=1<<m
    if type(k) is not int or not 1<=k<=min(m,4):raise ValueError('flat_rank out of range')
    context.limit('max_rank_work',comb(q-1,k)*q*(1<<k))
    spaces=set()
    for basis in combinations(range(1,q),k):
        H={0}
        for x in basis:H|={a^x for a in list(H)}
        if len(H)==1<<k:spaces.add(tuple(sorted(H)))
    flats=set()
    for H in spaces:
        unused=set(range(q))
        while unused:
            t=min(unused);S=tuple(sorted(t^x for x in H));flats.add(S);unused.difference_update(S)
    context.limit('max_codewords',len(flats))
    return sorted(flats)


def generate(inputs,params,context):
    check_params(params,('m','modulus','flat_rank','terms','exponents','coefficient_options','translations','seed_supports','max_seed_conflicts','dimension','coordinate_maps','allow_nonpermutations'))
    if inputs:raise ValueError('finite-field generation takes no inputs')
    m=params.get('m');mod=params.get('modulus');q=_field(m,mod);k=params.get('flat_rank',3)
    d=params.get('dimension',q);integer(d,'dimension');context.limit('max_dimension',d)
    if d<q:raise ValueError('dimension cannot be smaller than the field coordinate set')
    maps=params.get('coordinate_maps',[list(range(q))])
    if not isinstance(maps,list) or not maps or any(not isinstance(a,list) or len(a)!=q or any(type(x) is not int or not 0<=x<d for x in a) or len(set(a))!=q for a in maps):raise ValueError('coordinate maps must be injections into the ambient dimension')
    if len({tuple(a) for a in maps})!=len(maps):raise ValueError('duplicate coordinate maps')
    context.limit('max_trials',len(maps))
    terms=params.get('terms');exponents=params.get('exponents');options=params.get('coefficient_options')
    if terms is not None:
        if exponents is not None or options is not None or not isinstance(terms,list) or not terms:raise ValueError('provide terms or a coefficient grid')
        family=[terms];strict=True
    else:
        if not isinstance(exponents,list) or not exponents or not isinstance(options,list) or len(options)!=len(exponents) or any(not isinstance(v,list) or not v for v in options):raise ValueError('invalid coefficient grid')
        count=1
        for v in options:count*=len(v)
        context.limit('max_trials',count);family=[[list(z) for z in zip(coeff,exponents)] for coeff in product(*options)];strict=False
    allow_nonperm=params.get('allow_nonpermutations',False)
    if type(allow_nonperm) is not bool:raise ValueError('allow_nonpermutations must be boolean')
    maps_field=[];seenperm=set();invalid=0;permutation_count=0
    for terms in family:
        for term in terms:
            if not isinstance(term,list) or len(term)!=2 or any(type(x) is not int for x in term) or not(0<=term[0]<q) or not(0<=term[1]<=2*q):raise ValueError('invalid polynomial term')
        perm=[]
        for x in range(q):
            y=0
            for a,e in terms:y^=_mul(a,_power(x,e,m,mod),m,mod)
            perm.append(y)
        isperm=len(set(perm))==q
        if not isperm:
            invalid+=1
            if not allow_nonperm:
                if strict:raise ValueError('polynomial is not a permutation')
                continue
        if tuple(perm) not in seenperm:
            seenperm.add(tuple(perm));maps_field.append(perm)
            permutation_count+=int(isperm)
    if not maps_field:raise ValueError('no admissible polynomial map found')
    tr=params.get('translations',True)
    if type(tr) is not bool:raise ValueError('translations must be boolean')
    seed=params.get('seed_supports',[])
    if not isinstance(seed,list) or any(type(x) is not int or x<0 or x>>d or x.bit_count()!=(1<<k) for x in seed) or len(set(seed))!=len(seed):raise ValueError('invalid seed supports')
    maxconf=params.get('max_seed_conflicts',len(seed))
    if type(maxconf) is not int or maxconf<0:raise ValueError('max_seed_conflicts must be nonnegative')
    flats=affine_flats(m,k,context);context.limit('max_rank_work',len(flats)*(q if tr else 1)*len(maps_field)*len(maps));images=set();count=0
    collapsed=0
    for perm in maps_field:
        for S in flats:
            img=[perm[x] for x in S]
            if len(set(img))!=len(S):collapsed+=1;continue
            for t in (range(q) if tr else (0,)):
                for mapping in maps:
                    images.add(sum(1<<mapping[x^t] for x in img));count+=1
                    context.limit('max_codewords',len(images))
    context.limit('max_rank_work',len(images)*(len(seed)+1));seedset=set(seed);found=[]
    for mask in images:
        if mask not in seedset and sum((mask&old).bit_count()>(1<<k)//2 for old in seed)<=maxconf:found.append(mask)
    result={'backend':'finite_field_support_pool','parameters':params,'supports':sorted(found),'flats':len(flats),'images_examined':count,'unique_images':len(images),'polynomials_examined':len(family),'maps_considered':len(maps_field),'permutations':permutation_count,'nonpermutations':invalid,'collapsed_flats':collapsed,'scope':'algebraic support proposals, not a spherical code'}
    return Artifact('search_result',result,{},Verification('unverified','finite-field support generation; geometric assembly remains to be checked'))


def _lines(s,parities):
    pos=[j for j in range(s.bit_length()) if s>>j&1];assert len(pos)==8
    return [sum(1<<pos[j] for j in range(8) if z>>j&1) for z in range(256) if not(z&1) and z.bit_count()%2 in parities]


def _bad(s,v,t,w):
    ov=s&t
    return abs(ov.bit_count()-2*(ov&(v^w)).bit_count())>4


def exchange(inputs,params,context):
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix,vstack
    check_params(params,('seed_supports','exceptional_lines','parities','candidate_limit','seconds'))
    if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='finite_field_support_pool':raise ValueError('finite-field support pool required')
    pool=inputs[0].data['supports'];p=inputs[0].data['parameters'];q=1<<p['m'];d=p.get('dimension',q)
    if p.get('flat_rank',3)!=3:raise ValueError('signed octad exchange requires 8-element supports')
    if not isinstance(pool,list) or len(set(pool))!=len(pool) or any(type(s) is not int or s<0 or s>>d or s.bit_count()!=8 for s in pool):raise ValueError('invalid candidate supports')
    F=params.get('seed_supports',[]);E=params.get('exceptional_lines',[])
    if not isinstance(F,list) or len(set(F))!=len(F) or any(type(s) is not int or s<0 or s>>d or s.bit_count()!=8 for s in F):raise ValueError('invalid full seed supports')
    if not isinstance(E,list) or any(not isinstance(z,list) or len(z)!=2 or any(type(x) is not int for x in z) for z in E):raise ValueError('invalid exceptional lines')
    E=[tuple(z) for z in E]
    if len(set(E))!=len(E) or any(s<0 or s>>d or s.bit_count()!=8 or w&~s or w&(s&-s) for s,w in E):raise ValueError('invalid canonical exceptional line')
    if any(s in set(F) for s,w in E):raise ValueError('exceptional lines overlap full supports')
    parity=params.get('parities',[0,1])
    if not isinstance(parity,list) or not parity or len(set(parity))!=len(parity) or any(type(v) is not int or v not in (0,1) for v in parity):raise ValueError('parities must select 0 and/or 1')
    limit=params.get('candidate_limit',len(pool));integer(limit,'candidate_limit');C=pool[:limit]
    # A new line may reuse an existing support. Exclude exact source lines,
    # not the entire support: the other sign parity can still be searched.
    Fset=set(F);Eset=set(E);S=[]
    for s in C:S.extend((s,v) for v in _lines(s,parity) if not(s in Fset and v.bit_count()%2==0) and (s,v) not in Eset)
    M=len(S);context.limit('max_codewords',M);context.limit('max_pairs',M*(M-1)//2);context.limit('max_gram_entries',M*M)
    sec=params.get('seconds',1)
    if type(sec) not in (int,float) or not isfinite(sec) or sec<=0:raise ValueError('seconds must be positive')
    context.limit('solver_seconds',sec)
    # Preflight the actual support-level source scans and reuse their line
    # lists across all signs in a support. This avoids repeated family work.
    context.limit('max_rank_work',len(C)*(len(F)+len(E))+M*(M-1)//2)
    groups={}
    for s in C:
        oldlines=[]
        for t in F:
            if (s&t).bit_count()>4:oldlines.extend((t,w) for w in _lines(t,[0]))
        oldlines.extend((t,w) for t,w in E if (s&t).bit_count()>4)
        groups[s]=oldlines
    context.limit('max_rank_work',len(C)*(len(F)+len(E))+sum(sum(1 for t,v in S if t==s)*len(groups[s]) for s in C)+M*(M-1)//2)
    old={};a=[];b=[]
    for i,(s,v) in enumerate(S):
        for line in groups[s]:
            if _bad(s,v,*line):
                if line not in old:old[line]=len(old)
                a.append(i);b.append(old[line])
    I=[];J=[]
    for i in range(M):
        for j in range(i):
            if _bad(*S[i],*S[j]):I.append(j);J.append(i)
    N=len(old);context.limit('max_codewords',M+N);ne=len(a);e=len(I);context.limit('max_rank_work',max(1,ne+e)*4)
    rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)]
    cols=np.r_[np.c_[a,M+np.asarray(b,dtype=int)].ravel(),np.c_[I,J].ravel()]
    vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)]
    mat=coo_matrix((vals,(rows,cols)),shape=(ne+e,M+N)).tocsr();obj=np.r_[-np.ones(M),np.ones(N)]
    if not M:raise ValueError('empty candidate pool')
    r=milp(obj,integrality=np.ones(M+N),bounds=Bounds(0,1),constraints=LinearConstraint(mat,-np.inf,np.r_[np.zeros(ne),np.ones(e)]),options={'time_limit':float(sec),'presolve':True})
    chosen=[];removed=[]
    if r.x is not None:
        chosen=np.flatnonzero(r.x[:M]>.5).tolist();removed=np.flatnonzero(r.x[M:]>.5).tolist()
    oldkeys=[None]*N
    for key,j in old.items():oldkeys[j]=key
    oldremoved={oldkeys[j] for j in removed}
    for i in chosen:
        s,v=S[i]
        for t in F:
            if (s&t).bit_count()>4:
                for w in _lines(t,[0]):
                    if _bad(s,v,t,w) and (t,w) not in oldremoved:raise ValueError('unremoved old blocker')
        for t,w in E:
            if _bad(s,v,t,w) and (t,w) not in oldremoved:raise ValueError('unremoved exceptional blocker')
    if any(_bad(*S[i],*S[j]) for i,j in combinations(chosen,2)):raise ValueError('selected candidate conflict')
    gain=len(chosen)-len(removed)
    if gain<0:chosen=[];removed=[];gain=0
    data={'backend':'finite_field_signed_octad_exchange','candidate_lines':[list(S[i]) for i in chosen],
          'removed_lines':[list(oldkeys[j]) for j in removed], 'source_supports_sha256':digest(F),
          'source_exceptional_sha256':digest(E),'solver_status':str(r.message),
          'scope':'source validity and omitted source layers require independent certification'}
    return Artifact('exchange_proposal',data,{'gain_lines':gain,'gain_points_if_certified':2*gain},
                    Verification('unverified','candidate exchange checked; full source not independently certified'))


OPERATORS={'finite_field.support_pool':generate,'finite_field.octad_exchange':exchange}
INPUT_KINDS={'finite_field.support_pool':(),'finite_field.octad_exchange':('search_result',)}
OUTPUT_KINDS={'finite_field.support_pool':('search_result',),'finite_field.octad_exchange':('exchange_proposal',)}
