"""Joint numerical discovery and exact replay of several finite-group orbits.

The numerical phase searches orbits in a common subgroup fixed space and makes no proof claim.
Promotion reconstructs the complete exact rational matrix group and checks
all relative orbit products with integer arithmetic, including unequal seed
norms. This permits algebraic normalization without floating-point evidence.
"""
from collections import deque
from pathlib import Path
import math,time
from fractions import Fraction
from functools import reduce
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from ..core import Artifact,RunContext,Verification,check_params,integer,digest
from ..exact_gram import exact_integer_rank
from .exact_orbits import _generators


def _closure(raw,dimension,context):
    gens=_generators(raw,dimension,context);den=gens[0][1]
    if any(q!=den for m,q in gens):raise ValueError('generators must share a common denominator')
    if dimension*den*den>=2**62:raise ValueError('matrix products exceed safe exact integer range')
    gens=[np.asarray(m,dtype=np.int64) for m,q in gens];ident=np.eye(dimension,dtype=np.int64)*den
    groups=[ident];seen={tuple(ident.ravel())};queue=deque([ident])
    while queue:
        x=queue.popleft()
        for g in gens:
            y=g@x
            if np.any(y%den):raise ValueError('group closure leaves supplied common-denominator lattice')
            y=y//den;key=tuple(y.ravel())
            if key not in seen:
                context.limit('max_orbit_points',len(groups)+1);context.limit('max_gram_entries',(len(groups)+1)*dimension*dimension);seen.add(key);groups.append(y);queue.append(y)
    return np.asarray(groups,dtype=np.int64),den



def _fixed_integer_basis(group,den,d):
    """Exact rational nullspace of all stabilizer matrices minus identity."""
    rows=[]
    for M in group:
        delta=M-den*np.eye(d,dtype=np.int64)
        rows.extend([[Fraction(int(x)) for x in row] for row in delta if np.any(row)])
    pivots=[];r=0
    for j in range(d):
        pivot=next((i for i in range(r,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r];v=rows[r][j];rows[r]=[x/v for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][j]:
                q=rows[i][j];rows[i]=[x-q*y for x,y in zip(rows[i],rows[r])]
        pivots.append(j);r+=1
        if r==len(rows):break
    basis=[]
    for j in range(d):
        if j in pivots:continue
        v=[Fraction(0) for _ in range(d)];v[j]=Fraction(1)
        for i,p in enumerate(pivots):v[p]=-rows[i][j]
        denoms=[x.denominator for x in v];mult=math.lcm(*denoms);w=[int(x*mult) for x in v];g=reduce(math.gcd,(abs(x) for x in w));w=[x//g for x in w]
        if max(map(abs,w),default=0)>=2**50:raise ValueError('stabilizer fixed-space basis exceeds safe numerical range')
        basis.append(w)
    return np.asarray(basis,dtype=np.int64).T if basis else np.zeros((d,0),np.int64)


def _parameters(p,context):
    check_params(p,('generators','dimension','orbit_count','attempts','max_iterations','seed','seconds','scale','initial_seeds','stabilizer_generators'))
    d=integer(p.get('dimension'),'dimension');context.limit('max_dimension',d);s=integer(p.get('orbit_count'),'orbit_count');context.limit('max_nodes',s)
    matrices,den=_closure(p.get('generators'),d,context);context.limit('max_pairs',len(matrices)*s*(s+1)//2)
    supplied=p.get('stabilizer_generators',[])
    if not isinstance(supplied,list):raise ValueError('stabilizer_generators must be a list')
    if supplied:
        subgroup,sd=_closure(supplied,d,context)
        if sd!=den:raise ValueError('stabilizer denominator must equal group denominator')
        group_set={tuple(M.ravel()) for M in matrices}
        if any(tuple(M.ravel()) not in group_set for M in subgroup):raise ValueError('stabilizer is not a subgroup of the group')
    else:subgroup=np.asarray([np.eye(d,dtype=np.int64)*den])
    sub_set={tuple(M.ravel()) for M in subgroup};exclude=np.array([tuple(M.ravel()) in sub_set for M in matrices],bool)
    basis=_fixed_integer_basis(subgroup,den,d)
    if basis.shape[1]==0:raise ValueError('stabilizer has no nonzero fixed vectors')
    return matrices,den,d,s,exclude,basis


def search_group_seeds(p,context=None):
    context=context or RunContext(Path('.'));M,den,d,s,exclude,basis=_parameters(p,context);attempts=integer(p.get('attempts',5),'attempts');context.limit('max_trials',attempts)
    maxiter=integer(p.get('max_iterations',200),'max_iterations');context.limit('max_rank_work',len(M)*d*d*s*maxiter)
    scale=integer(p.get('scale',10000),'scale');context.limit('max_rational_bits',scale.bit_length())
    if scale>100000000:raise ValueError('scale exceeds safe numerical rationalization range')
    seed=integer(p.get('seed',0),'seed',0);seconds=p.get('seconds',10)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);initial=p.get('initial_seeds')
    if initial is not None:
        if not isinstance(initial,list) or len(initial)!=s or any(not isinstance(row,list) or len(row)!=d or any(type(v) not in (int,float) or not math.isfinite(v) for v in row) for row in initial):raise ValueError('invalid numerical initial seeds')
        if any(sum(float(v)**2 for v in row)==0 for row in initial):raise ValueError('initial seeds must be nonzero')
    U,_=np.linalg.qr(basis.astype(float));r=U.shape[1];A=np.einsum('ia,gij,jb->gab',U,M.astype(float)/den,U,optimize=True);rng=np.random.default_rng(seed);deadline=time.monotonic()+seconds;best=math.inf;bestX=None;finished=0
    class Cutoff(Exception):pass
    def evaluate(v,beta,gradient=True):
        nonlocal best,bestX
        if time.monotonic()>deadline:raise Cutoff()
        Y=v.reshape(s,r);norm=np.linalg.norm(Y,axis=1)
        if np.any(norm<1e-14):raise ValueError('optimizer produced a zero seed')
        X=Y/norm[:,None];forward=np.einsum('gij,sj->sgi',A,X,optimize=True);back=np.einsum('gji,sj->sgi',A,X,optimize=True);vals=[];meta=[]
        for i in range(s):
            for j in range(i,s):
                z=forward[j]@X[i]
                if i==j:z=z[~exclude]
                if len(z):vals.append(z);meta.append((i,j))
        if not vals:raise ValueError('empty generic orbit objective')
        flat=np.concatenate(vals);maximum=float(flat.max())
        if maximum<best:best=maximum;bestX=X.copy()
        weight=np.exp(beta*flat-logsumexp(beta*flat));value=logsumexp(beta*flat)/beta;grad=np.zeros_like(X);at=0
        for z,(i,j) in zip(vals,meta):
            w=weight[at:at+len(z)];at+=len(z)
            if i==j:
                f=forward[i,~exclude];b=back[i,~exclude];grad[i]+=w@(f+b)
            else:grad[i]+=w@forward[j];grad[j]+=w@back[i]
        for i in range(s):grad[i]=(grad[i]-X[i]*(X[i]@grad[i]))/norm[i]
        return value,grad.ravel()
    for trial in range(attempts):
        if time.monotonic()>=deadline:break
        Y=np.asarray(initial,float)@U if trial==0 and initial is not None else rng.normal(size=(s,r))
        if np.any(np.linalg.norm(Y,axis=1)<1e-14):raise ValueError('initial seed projects to zero in stabilizer fixed space')
        try:
            for beta in (12.,35.,90.,180.):
                result=minimize(lambda v:evaluate(v,beta),Y.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-10,'gtol':1e-7});Y=result.x.reshape(s,r)
            finished+=1
        except Cutoff:break
    if bestX is None:return {'parameters':p,'group_order':len(M),'finished_attempts':finished,'status':'no numerical iterate','seed_vectors':[]}
    vectors=bestX@U.T;coeff=np.linalg.lstsq(basis.astype(float),vectors.T,rcond=None)[0].T;scaled=np.rint(coeff*scale)
    if not np.all(np.isfinite(scaled)) or float(np.max(np.abs(scaled),initial=0))*int(np.max(np.abs(basis),initial=0))*r>=2**62:raise ValueError('rationalized seeds exceed safe integer range')
    seeds=scaled.astype(np.int64)@basis.T
    if int(np.max(np.abs(seeds),initial=0))*d*den>=2**62:raise ValueError('stabilizer checks exceed safe integer range')
    if not np.all(np.einsum('gij,sj->gsi',M[exclude],seeds)==den*seeds[None,:,:]):raise ValueError('rounded seeds do not lie in the fixed space')
    return {'parameters':p,'group_order':len(M),'finished_attempts':finished,'status':'numerical proposal only','best_max_correlation':best,
            'seed_vectors':seeds.tolist(),'stabilizer_order':int(exclude.sum()),'orbit_cardinality_upper_forecast':s*len(M)//int(exclude.sum()),
            'scope':'Numerical orbit search in a common stabilizer fixed space. Additional stabilizers may reduce cardinality; promote only after exact verification.'}


def group_seed_search(inputs,params,context):
    if inputs:raise ValueError('group seed search takes no inputs')
    result=search_group_seeds(params,context)
    return Artifact('search_result',{'backend':'group_seed_search','format':'group_seed_search_v1','result':result},{},Verification('unverified','numerical multi-orbit proposal; exact promotion required',{}))


def verify_search_artifact(artifact,context=None):
    context=context or RunContext(Path('.'));data=artifact.data
    if artifact.kind!='search_result' or data.get('backend')!='group_seed_search' or data.get('format')!='group_seed_search_v1':raise ValueError('invalid group search receipt')
    r=data.get('result');p=r.get('parameters') if isinstance(r,dict) else None
    if not isinstance(p,dict):raise ValueError('missing parameters')
    M,den,d,s,exclude,basis=_parameters(p,context);seeds=r.get('seed_vectors')
    if not isinstance(seeds,list) or (seeds and (len(seeds)!=s or any(not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row) for row in seeds))):raise ValueError('invalid proposal seeds')
    for row in seeds:
        for value in row:context.limit('max_rational_bits',abs(value).bit_length())
    return Verification('unverified','exact group generators replayed; numerical proposal not certified',{'facts':{},'group_order':len(M)})


def verify_multi_artifact(artifact,context=None):
    context=context or RunContext(Path('.'));data=artifact.data
    if artifact.kind!='implicit_spherical_code' or not isinstance(data,dict) or set(data)!={'backend','generators','seed_vectors'} or data.get('backend')!='exact_multi_group_orbits':raise ValueError('invalid exact multi-orbit payload')
    seeds=data['seed_vectors']
    if not isinstance(seeds,list) or not seeds or not isinstance(seeds[0],list) or not seeds[0]:raise ValueError('missing seeds')
    d=len(seeds[0]);context.limit('max_dimension',d);context.limit('max_nodes',len(seeds))
    if any(not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row) for row in seeds):raise ValueError('invalid exact seeds')
    M,den=_closure(data['generators'],d,context);s=len(seeds);context.limit('max_pairs',len(M)*s*(s+1)//2)
    maxseed=max(abs(x) for row in seeds for x in row);context.limit('max_rational_bits',maxseed.bit_length())
    if not maxseed or d*den*maxseed>=2**62:raise ValueError('seed multiplication out of range')
    orbits=[];norms=[]
    for row in seeds:
        p=np.asarray(row,dtype=np.int64);norm=den*den*sum(int(x)*int(x) for x in row)
        if not norm:raise ValueError('zero seed')
        V=M@p;unique={tuple(map(int,x)) for x in V};context.limit('max_materialized_points',sum(map(len,orbits))+len(unique));orbits.append(unique);norms.append(norm)
    # Transitivity of the exactly closed orthogonal group reduces all pairs
    # to one identity representative against every relative group image.
    for i in range(s):
        seed=tuple(den*x for x in seeds[i]);ni=norms[i]
        for j in range(i,s):
            nj=norms[j]
            for v in orbits[j]:
                if i==j and v==seed:continue
                dot=sum(int(x)*int(y) for x,y in zip(seed,v))
                if dot>0 and 4*dot*dot>ni*nj:raise ValueError('exact cross-orbit separation fails')
    rows=[]
    for orbit in orbits:
        rows.extend([list(v) for v in orbit])
        if len(rows)>=max(64,4*d):break
    rank,_=exact_integer_rank(rows,context)
    if rank<d:
        rows=[list(v) for orbit in orbits for v in orbit];rank,_=exact_integer_rank(rows,context)
    facts={'dimension':rank,'rank_exact':rank,'ambient_dimension':d,'cardinality':sum(map(len,orbits)),
           'max_inner_product_upper_bound':'1/2'}
    return Verification('exact','closed rational group and exact relative-orbit product inequalities',
                        {'facts':facts,'group_order':len(M),'orbit_sizes':list(map(len,orbits)),'data_sha256':digest(data)})


def certify_group_seeds(inputs,params,context):
    check_params(params,())
    if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='group_seed_search':raise ValueError('expected a group seed search result')
    verify_search_artifact(inputs[0],context);r=inputs[0].data['result'];seeds=r['seed_vectors']
    if not seeds:raise ValueError('search returned no seeds')
    a=Artifact('implicit_spherical_code',{'backend':'exact_multi_group_orbits','generators':r['parameters']['generators'],'seed_vectors':seeds},{},Verification())
    a.verification=verify_multi_artifact(a,context);a.facts=a.verification.details['facts'];return a

OPERATORS={'search.group_orbit_seeds':group_seed_search,'groups.certify_multi_orbit':certify_group_seeds}
INPUT_KINDS={'search.group_orbit_seeds':(),'groups.certify_multi_orbit':('search_result',)}
OUTPUT_KINDS={'search.group_orbit_seeds':('search_result',),'groups.certify_multi_orbit':('implicit_spherical_code',)}
