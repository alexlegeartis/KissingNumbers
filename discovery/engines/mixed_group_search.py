"""Numerically search group-orbit unions with an individual stabilizer per seed.

Results remain proposals; exact independent replay reconstructs and checks all
orbits. Different subgroups must lie in one exact rational orthogonal group.
"""
import math,time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from ..core import Artifact,RunContext,Verification,integer,check_params
from .group_seed_search import _closure,_fixed_integer_basis


def _setup(params,context):
    allowed=('generators','dimension','orbit_specs','attempts','max_iterations','seed','seconds','scale','initial_seeds','keep_initial')
    check_params(params,allowed);d=integer(params.get('dimension'),'dimension');context.limit('max_dimension',d)
    specs=params.get('orbit_specs')
    if not isinstance(specs,list) or not specs:raise ValueError('orbit_specs must be a nonempty list')
    s=len(specs);context.limit('max_nodes',s)
    M,den=_closure(params.get('generators'),d,context);context.limit('max_pairs',len(M)*s*(s+1)//2)
    keys={tuple(g.ravel()) for g in M};bases=[];omit=[];orders=[]
    for spec in specs:
        check_params(spec,('stabilizer_generators',));g=spec.get('stabilizer_generators',[])
        if not isinstance(g,list):raise ValueError('stabilizer_generators must be a list')
        if g:
            H,hd=_closure(g,d,context)
            if hd!=den:raise ValueError('stabilizer denominator must equal group denominator')
            if any(tuple(h.ravel()) not in keys for h in H):raise ValueError('stabilizer is not a subgroup')
        else:H=np.array([den*np.eye(d,dtype=np.int64)])
        hs={tuple(h.ravel()) for h in H};basis=_fixed_integer_basis(H,den,d)
        if basis.shape[1]==0:raise ValueError('stabilizer has no fixed vectors')
        bases.append(basis);orders.append(len(H));omit.append(np.array([tuple(h.ravel()) not in hs for h in M],bool))
    return M,den,d,s,bases,omit,orders


def search_mixed(params,context=None):
    context=context or RunContext(__import__('pathlib').Path('.'));M,den,d,s,bases,include,orders=_setup(params,context)
    attempts=integer(params.get('attempts',5),'attempts');maxiter=integer(params.get('max_iterations',200),'max_iterations');scale=integer(params.get('scale',10000),'scale')
    context.limit('max_trials',attempts);context.limit('max_rank_work',len(M)*d*d*s*maxiter);context.limit('max_rational_bits',scale.bit_length())
    if scale>100000000:raise ValueError('scale exceeds safe range')
    seconds=params.get('seconds',10)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);seed=integer(params.get('seed',0),'seed',0);rng=np.random.default_rng(seed)
    initial=params.get('initial_seeds')
    if initial is not None and (not isinstance(initial,list) or len(initial)!=s or any(not isinstance(v,list) or len(v)!=d or any(type(x) not in (int,float) or not math.isfinite(x) for x in v) for v in initial)):raise ValueError('invalid initial seeds')
    keep=params.get('keep_initial',False)
    if type(keep) is not bool:raise ValueError('keep_initial must be boolean')
    if keep and (initial is None or any(type(x) is not int for row in initial for x in row)):raise ValueError('keep_initial requires integer initial seeds')
    if keep:
        for i,row in enumerate(initial):
            if max(map(abs,row),default=0)*d*den>=2**62:raise ValueError('initial fixed-space check exceeds safe range')
            v=np.array(row,dtype=np.int64)
            if not np.any(v) or not np.all(np.einsum('gij,j->gi',M[~include[i]],v)==den*v):raise ValueError('initial seed is not fixed by its stabilizer')
    U=[np.linalg.qr(b.astype(float))[0] for b in bases];offset=np.cumsum([0]+[u.shape[1] for u in U]);A=M.astype(float)/den;best=math.inf;bestX=None;deadline=time.monotonic()+seconds;finished=0
    class Cutoff(Exception):pass
    def evaluate(y,beta):
        nonlocal best,bestX
        if time.monotonic()>deadline:raise Cutoff()
        X=[];norm=[]
        for i,u in enumerate(U):
            q=y[offset[i]:offset[i+1]];n=np.linalg.norm(q)
            if n<1e-14:raise ValueError('zero seed in fixed space')
            norm.append(n);X.append(u@(q/n))
        X=np.asarray(X);forward=np.einsum('gij,sj->sgi',A,X,optimize=True);back=np.einsum('gji,sj->sgi',A,X,optimize=True);vals=[];meta=[]
        for i in range(s):
            for j in range(i,s):
                z=forward[j]@X[i]
                if i==j:z=z[include[i]]
                if len(z):vals.append(z);meta.append((i,j))
        if not vals:raise ValueError('empty orbit objective')
        flat=np.concatenate(vals);mx=float(flat.max())
        if mx<best:best=mx;bestX=X.copy()
        lse=logsumexp(beta*flat);weight=np.exp(beta*flat-lse);grad=np.zeros_like(X);at=0
        for z,(i,j) in zip(vals,meta):
            w=weight[at:at+len(z)];at+=len(z)
            if i==j:grad[i]+=w@(forward[i,include[i]]+back[i,include[i]])
            else:grad[i]+=w@forward[j];grad[j]+=w@back[i]
        gy=np.empty(offset[-1])
        for i,u in enumerate(U):gy[offset[i]:offset[i+1]]=u.T@(grad[i]-X[i]*(X[i]@grad[i]))/norm[i]
        return lse/beta,gy
    for trial in range(attempts):
        if time.monotonic()>=deadline:break
        if trial==0 and initial is not None:y=np.concatenate([u.T@np.asarray(v,float) for u,v in zip(U,initial)])
        else:y=rng.normal(size=offset[-1])
        if any(np.linalg.norm(y[offset[i]:offset[i+1]])<1e-14 for i in range(s)):raise ValueError('initial seed projects to zero')
        try:
            for beta in (12.,35.,90.,180.):
                res=minimize(lambda v:evaluate(v,beta),y,jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-10,'gtol':1e-7});y=res.x
            finished+=1
        except Cutoff:break
    if bestX is None:return {'parameters':{'generators':params['generators'],'dimension':d,'orbit_count':s},'seed_vectors':[],'status':'no numerical iterate'}
    seeds=[]
    for i,b in enumerate(bases):
        coeff=np.linalg.lstsq(b.astype(float),bestX[i],rcond=None)[0];q=np.rint(coeff*scale)
        if not np.all(np.isfinite(q)) or float(np.max(np.abs(q),initial=0))*int(np.max(np.abs(b),initial=0))*len(q)>=2**62:raise ValueError('rationalized seeds exceed safe range')
        v=q.astype(np.int64)@b.T
        if int(np.max(np.abs(v),initial=0))*d*den>=2**62:raise ValueError('stabilizer check exceeds safe range')
        if not np.all(np.einsum('gij,j->gi',M[~include[i]],v)==den*v):raise ValueError('rationalized seed leaves its fixed space')
        seeds.append(v.tolist())
    # Optional exact initial values are retained if a numerical search did not
    # establish strict slack. Promotion still checks all orbit products.
    if keep and best>=.5-1e-8:seeds=initial
    p={'generators':params['generators'],'dimension':d,'orbit_count':s}
    return {'parameters':p,'seed_vectors':seeds,'group_order':len(M),'stabilizer_orders':orders,
            'orbit_cardinality_upper_forecast':sum(len(M)//n for n in orders),
            'best_max_correlation':best,'finished_attempts':finished,'status':'numerical proposal only',
            'scope':'Per-orbit stabilizers; orbit sizes and all distances require independent exact promotion.'}


def mixed_search(inputs,params,context):
    if inputs:raise ValueError('mixed stabilizer search takes no inputs')
    r=search_mixed(params,context)
    return Artifact('search_result',{'backend':'group_seed_search','format':'group_seed_search_v1','result':r},{},Verification('unverified','mixed-stabilizer proposal; exact promotion required',{}))

OPERATORS={'search.mixed_stabilizer_orbits':mixed_search}
INPUT_KINDS={'search.mixed_stabilizer_orbits':()}
OUTPUT_KINDS={'search.mixed_stabilizer_orbits':('search_result',)}
