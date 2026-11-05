"""Jointly optimize orbits of different groups, including singleton points.

Unlike a single-group orbit union, this permits deliberate symmetry breaking.
All numerical output is provisional; exact promotion materializes every orbit
and checks every pair with integer arithmetic, including unequal norms.
"""
from pathlib import Path
import math,time
import numpy as np
from scipy.optimize import minimize
from ..core import Artifact,RunContext,Verification,integer,check_params,digest
from ..exact_gram import exact_integer_rank
from .group_seed_search import _closure,_fixed_integer_basis
from .cavity_repacking import _smooth_value_gradient


def _group(g,d,context,identity_den=1):
    if not isinstance(g,list):raise ValueError('generators must be a list')
    if not g:return np.array([np.eye(d,dtype=np.int64)*identity_den]),identity_den
    return _closure(g,d,context)


def _local(spec,d,context):
    check_params(spec,('generators','stabilizer_generators'))
    M,den=_group(spec.get('generators',[]),d,context);H,hd=_group(spec.get('stabilizer_generators',[]),d,context,den)
    if hd!=den:raise ValueError('stabilizer denominator differs from orbit group')
    keys={tuple(g.ravel()):i for i,g in enumerate(M)}
    if any(tuple(h.ravel()) not in keys for h in H):raise ValueError('stabilizer is not a subgroup of orbit group')
    owner=set();reps=[]
    for i,g in enumerate(M):
        if i in owner:continue
        reps.append(g)
        for h in H:
            q=g@h
            if np.any(q%den):raise ValueError('orbit group leaves its exact rational lattice')
            key=tuple((q//den).ravel())
            if key not in keys or keys[key] in owner:raise ValueError('invalid orbit coset decomposition')
            owner.add(keys[key])
    if len(owner)!=len(M):raise ValueError('incomplete orbit coset decomposition')
    basis=_fixed_integer_basis(H,den,d)
    if basis.shape[1]==0:raise ValueError('stabilizer has no fixed vectors')
    return M,den,H,basis,np.array(reps,dtype=np.int64)


def _params(p,context):
    check_params(p,('dimension','orbit_specs','attempts','max_iterations','seed','seconds','scale','initial_seeds','keep_initial','betas','metric'))
    d=integer(p.get('dimension'),'dimension');context.limit('max_dimension',d);specs=p.get('orbit_specs')
    if not isinstance(specs,list) or not specs:raise ValueError('orbit_specs must be nonempty')
    context.limit('max_nodes',len(specs));parts=[_local(spec,d,context) for spec in specs];N=sum(len(t[4]) for t in parts)
    context.limit('max_materialized_points',N);context.limit('max_pairs',N*(N-1)//2);context.limit('max_gram_entries',N*N)
    metric=p.get('metric',[1]*d)
    if not isinstance(metric,list) or len(metric)!=d or any(type(w) is not int or w<=0 for w in metric):raise ValueError('metric must be positive integer diagonal')
    if max(metric)*d*max(t[1] for t in parts)**2>=2**62:raise ValueError('metric group check exceeds safe range')
    W=np.diag(np.asarray(metric,dtype=np.int64))
    for M,den,H,basis,reps in parts:
        if not np.all(np.einsum('gji,jk,gkl->gil',M,W,M)==den*den*W):raise ValueError('group does not preserve the metric')
    return d,specs,parts,N,metric


def search_partial(p,context):
    d,specs,parts,N,metric=_params(p,context);s=len(specs);attempts=integer(p.get('attempts',5),'attempts');maxiter=integer(p.get('max_iterations',100),'max_iterations');scale=integer(p.get('scale',10000),'scale');context.limit('max_trials',attempts)
    context.limit('max_rank_work',N*N*maxiter);context.limit('max_rational_bits',scale.bit_length())
    if scale>100000000:raise ValueError('scale exceeds safe range')
    betas=p.get('betas',[12.,35.,90.,180.])
    if not isinstance(betas,list) or not 1<=len(betas)<=16 or any(type(b) not in (int,float) or not math.isfinite(b) or b<=0 or b>100000 for b in betas):raise ValueError('invalid beta schedule')
    seconds=p.get('seconds',10)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);rng=np.random.default_rng(integer(p.get('seed',0),'seed',0));initial=p.get('initial_seeds');keep=p.get('keep_initial',False)
    if type(keep) is not bool:raise ValueError('keep_initial must be boolean')
    if initial is not None and (not isinstance(initial,list) or len(initial)!=s or any(not isinstance(v,list) or len(v)!=d or any(type(x) not in (int,float) or not math.isfinite(x) for x in v) for v in initial)):raise ValueError('invalid initial seeds')
    if keep and (initial is None or any(type(x) is not int for row in initial for x in row)):raise ValueError('keep_initial requires integer seeds')
    if keep:
        for row,(M,den,H,basis,reps) in zip(initial,parts):
            if max(map(abs,row),default=0)*d*den>=2**62:raise ValueError('initial seed is too large')
            v=np.array(row,dtype=np.int64)
            if not np.any(v) or not np.all(np.einsum('gij,j->gi',H,v)==den*v):raise ValueError('initial seed is not fixed by its stabilizer')
    sqrtw=np.sqrt(np.asarray(metric,float));U=[np.linalg.qr(sqrtw[:,None]*t[3].astype(float))[0] for t in parts];off=np.cumsum([0]+[u.shape[1] for u in U]);trans=[sqrtw[None,:,None]*t[4].astype(float)/t[1]/sqrtw[None,None,:] for t in parts];span=np.cumsum([0]+[len(x) for x in trans]);deadline=time.monotonic()+seconds;best=math.inf;bestX=None;finished=0
    class Cutoff(Exception):pass
    def evaluate(y,beta):
        nonlocal best,bestX
        if time.monotonic()>deadline:raise Cutoff()
        X=[];norm=[]
        for i,u in enumerate(U):
            z=y[off[i]:off[i+1]];q=np.linalg.norm(z)
            if q<1e-14:raise ValueError('seed projects to zero')
            norm.append(q);X.append(u@(z/q))
        X=np.array(X);V=np.concatenate([np.einsum('gij,j->gi',t,x) for t,x in zip(trans,X)])
        value,gradV,mx=_smooth_value_gradient(V,np.empty((0,d)),beta)
        if mx<best:best=mx;bestX=X.copy()
        gy=np.empty(off[-1])
        for i,u in enumerate(U):
            grad=np.einsum('gji,gj->i',trans[i],gradV[span[i]:span[i+1]])
            gy[off[i]:off[i+1]]=u.T@(grad-X[i]*(X[i]@grad))/norm[i]
        return value,gy
    for trial in range(attempts):
        if time.monotonic()>=deadline:break
        y=np.concatenate([u.T@(sqrtw*np.array(v,float)) for u,v in zip(U,initial)]) if trial==0 and initial is not None else rng.normal(size=off[-1])
        try:
            for beta in betas:
                r=minimize(lambda x:evaluate(x,beta),y,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-10,'gtol':1e-7});y=r.x
            finished+=1
        except Cutoff:break
    if bestX is None:return {'parameters':p,'seed_vectors':[],'status':'no numerical iterate'}
    seeds=[]
    for i,(_,den,H,basis,_) in enumerate(parts):
        q=np.rint(np.linalg.lstsq(basis.astype(float),bestX[i]/sqrtw,rcond=None)[0]*scale)
        if not np.all(np.isfinite(q)) or float(np.max(np.abs(q),initial=0))*int(np.max(np.abs(basis),initial=0))*len(q)>=2**62:raise ValueError('rationalized seed exceeds safe range')
        v=q.astype(np.int64)@basis.T
        if int(np.max(np.abs(v),initial=0))*d*den>=2**62 or not np.all(np.einsum('gij,j->gi',H,v)==den*v):raise ValueError('rationalized seed is not fixed')
        seeds.append(v.tolist())
    if keep and best>=.5-1e-8:seeds=initial
    return {'parameters':p,'seed_vectors':seeds,'best_max_correlation':best,'finished_attempts':finished,
            'orbit_cardinality_upper_forecast':N,'status':'numerical proposal only'}


def partial_search(inputs,params,context):
    if inputs:raise ValueError('partial-group search takes no inputs')
    return Artifact('search_result',{'backend':'partial_group_search','result':search_partial(params,context)},{},Verification('unverified','partial-orbit proposal; exact promotion required',{}))


def verify_search_artifact(a,context=None):
    context=context or RunContext(Path('.'));r=a.data.get('result',{});p=r.get('parameters');d,specs,parts,N,metric=_params(p,context);seeds=r.get('seed_vectors')
    if not isinstance(seeds,list) or (seeds and (len(seeds)!=len(specs) or any(not isinstance(row,list) or len(row)!=d or any(type(v) is not int for v in row) for row in seeds))):raise ValueError('invalid partial-orbit seeds')
    return Verification('unverified','partial orbit groups replayed; separation not certified',{'facts':{}})


def verify_exact(a,context=None):
    context=context or RunContext(Path('.'));data=a.data
    if a.kind!='implicit_spherical_code' or not isinstance(data,dict) or set(data)!={'backend','groups','seed_vectors','metric'} or data['backend']!='exact_partial_group_orbits':raise ValueError('invalid exact partial-orbit payload')
    seeds=data['seed_vectors'];groups=data['groups']
    if not isinstance(seeds,list) or not seeds or not isinstance(groups,list) or len(groups)!=len(seeds) or not isinstance(seeds[0],list) or not seeds[0]:raise ValueError('missing orbit data')
    d=len(seeds[0]);context.limit('max_dimension',d);context.limit('max_nodes',len(seeds));parts=[];count=0
    metric=data['metric']
    if not isinstance(metric,list) or len(metric)!=d or any(type(w) is not int or w<=0 for w in metric):raise ValueError('invalid exact metric')
    W=np.diag(np.asarray(metric,dtype=np.int64))
    for gens,row in zip(groups,seeds):
        if not isinstance(row,list) or len(row)!=d or any(type(v) is not int for v in row):raise ValueError('invalid integer seed')
        M,den=_group(gens,d,context)
        if max(metric)*d*den*den>=2**62 or not np.all(np.einsum('gji,jk,gkl->gil',M,W,M)==den*den*W):raise ValueError('group does not preserve the exact metric')
        mx=max(map(abs,row),default=0);context.limit('max_rational_bits',mx.bit_length())
        if not mx or mx*d*den>=2**62:raise ValueError('seed multiplication out of range')
        V=np.unique(np.einsum('gij,j->gi',M,np.asarray(row,dtype=np.int64)),axis=0);norm=den*den*sum(w*v*v for w,v in zip(metric,row));count+=len(V);context.limit('max_materialized_points',count);parts.append((V,norm))
    context.limit('max_pairs',count*(count-1)//2);batch=max(1,min(128,context.budget.get('max_gram_entries',2000000)//max(1,count)));context.limit('max_gram_entries',batch*count)
    for i,(A,na) in enumerate(parts):
        ma=int(np.max(np.abs(A),initial=0))
        for j in range(i,len(parts)):
            B,nb=parts[j];mb=int(np.max(np.abs(B),initial=0))
            if d*ma*mb*max(metric)>=2**62:raise ValueError('exact pair products exceed safe range')
            # For integer dot products, dot <= floor(sqrt(na*nb)/2) is
            # equivalent to 4*dot^2 <= na*nb when dot is positive.
            threshold=math.isqrt(na*nb)//2
            for start in range(0,len(A),batch):
                D=(A[start:start+batch]*np.asarray(metric,dtype=np.int64))@B.T
                if i==j:
                    for k in range(len(D)):D[k,start+k]=0
                if np.any(D>threshold):raise ValueError('exact partial-orbit separation fails')
    rows=[]
    for V,_ in parts:
        rows.extend(V.tolist())
        if len(rows)>=max(64,4*d):break
    rank,_=exact_integer_rank(rows,context)
    if rank<d:rank,_=exact_integer_rank([v.tolist() for V,_ in parts for v in V],context)
    facts={'dimension':rank,'rank_exact':rank,'ambient_dimension':d,'cardinality':count,'max_inner_product_upper_bound':'1/2'}
    return Verification('exact','materialized rational subgroup orbits and integer pair checks',{'facts':facts,'orbit_sizes':[len(x) for x,_ in parts],'data_sha256':digest(data)})


def certify_partial(inputs,params,context):
    check_params(params,())
    if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='partial_group_search':raise ValueError('expected partial-group proposal')
    verify_search_artifact(inputs[0],context);r=inputs[0].data['result'];seeds=r['seed_vectors']
    if not seeds:raise ValueError('search returned no seeds')
    data={'backend':'exact_partial_group_orbits','groups':[v.get('generators',[]) for v in r['parameters']['orbit_specs']],'seed_vectors':seeds,'metric':r['parameters'].get('metric',[1]*r['parameters']['dimension'])};a=Artifact('implicit_spherical_code',data,{},Verification());a.verification=verify_exact(a,context);a.facts=a.verification.details['facts'];return a

OPERATORS={'search.partial_group_assemblies':partial_search,'groups.certify_partial_assemblies':certify_partial}
INPUT_KINDS={'search.partial_group_assemblies':(),'groups.certify_partial_assemblies':('search_result',)}
OUTPUT_KINDS={'search.partial_group_assemblies':('search_result',),'groups.certify_partial_assemblies':('implicit_spherical_code',)}
