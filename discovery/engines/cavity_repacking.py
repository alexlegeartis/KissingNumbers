"""Jointly replace and move a bounded neighborhood around a spherical code.

The search optimizes all moving points against each other and against fixed
points, with an optional positive diagonal metric. Its numerical output is
never a bound until all fixed and moved points are exactly checked.
"""
from pathlib import Path
import math,time
import numpy as np
from scipy.optimize import minimize
from ..core import Artifact,RunContext,Verification,integer,check_params,digest
from ..exact_gram import exact_integer_rank


def _smooth_value_gradient(points,fixed,beta):
    """Log-sum-exp over each unordered moving pair and each moving/fixed pair.

    The symmetric moving Gram is weighted by one half and overwritten in
    place. This avoids explicit triangular indices and a second dense matrix.
    """
    D=points@fixed.T;R=points@points.T;np.fill_diagonal(R,-np.inf)
    mx=max(float(np.max(D,initial=-np.inf)),float(np.max(R,initial=-np.inf)))
    if not math.isfinite(mx):raise ValueError('need at least two total points')
    D-=mx;D*=beta;np.exp(D,out=D)
    R-=mx;R*=beta;np.exp(R,out=R);R*=.5
    total=float(D.sum()+R.sum())
    grad=(D@fixed+2*(R@points))/total
    return mx+math.log(total)/beta,grad,mx

def _source(p,context):
    check_params(p,('source','metric','remove_indices','initial_replacements','attempts','max_iterations','seed','seconds','scale','betas'))
    source=p.get('source');rm=p.get('remove_indices');initial=p.get('initial_replacements')
    if not isinstance(source,list) or not source or not isinstance(source[0],list) or not source[0]:raise ValueError('source must be a nonempty integer matrix')
    n=len(source);d=len(source[0]);context.limit('max_dimension',d);context.limit('max_materialized_points',n)
    if any(not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row) for row in source):raise ValueError('source coordinates must be integers')
    if not isinstance(rm,list) or any(type(i) is not int or i<0 or i>=n for i in rm) or len(set(rm))!=len(rm):raise ValueError('invalid removal indices')
    if not isinstance(initial,list) or not initial or any(not isinstance(v,list) or len(v)!=d or any(type(x) not in (int,float) or not math.isfinite(x) for x in v) for v in initial):raise ValueError('invalid replacement seeds')
    metric=p.get('metric',[1]*d)
    if not isinstance(metric,list) or len(metric)!=d or any(type(w) is not int or w<=0 for w in metric):raise ValueError('invalid positive diagonal metric')
    context.limit('max_nodes',len(initial));total=n-len(rm)+len(initial);context.limit('max_materialized_points',total);context.limit('max_pairs',total*(total-1)//2)
    context.limit('max_gram_entries',len(initial)*(n-len(rm))+len(initial)**2)
    return np.asarray(source,dtype=float),np.asarray(metric,dtype=float),set(rm),np.asarray(initial,dtype=float),d,total


def search_cavity(inputs,p,context):
    if inputs:raise ValueError('cavity search takes no inputs')
    source,metric,removed,initial,d,total=_source(p,context);m=len(initial);fixed=source[[i not in removed for i in range(len(source))]];w=np.sqrt(metric);fixed=fixed*w;fn=np.linalg.norm(fixed,axis=1)
    if np.any(fn<1e-15):raise ValueError('source contains zero vector')
    fixed/=fn[:,None];attempts=integer(p.get('attempts',1),'attempts');maxiter=integer(p.get('max_iterations',200),'max_iterations');scale=integer(p.get('scale',10000),'scale');context.limit('max_trials',attempts);context.limit('max_rational_bits',scale.bit_length())
    if scale>100000000:raise ValueError('scale exceeds safe range')
    context.limit('max_rank_work',(m*len(fixed)+m*m)*maxiter)
    betas=p.get('betas',[100.,300.,1000.,3000.])
    if not isinstance(betas,list) or not 1<=len(betas)<=16 or any(type(b) not in (int,float) or not math.isfinite(b) or b<=0 or b>100000 for b in betas):raise ValueError('invalid beta schedule')
    seconds=p.get('seconds',10)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
    context.limit('solver_seconds',seconds);rng=np.random.default_rng(integer(p.get('seed',0),'seed',0));best=math.inf;bestX=None;finished=0;deadline=time.monotonic()+seconds
    class Cutoff(Exception):pass
    def evaluate(y,beta):
        nonlocal best,bestX
        if time.monotonic()>deadline:raise Cutoff()
        Y=y.reshape(m,d);norm=np.linalg.norm(Y,axis=1)
        if np.any(norm<1e-14):raise ValueError('zero replacement seed')
        Z=Y/norm[:,None];value,grad,mx=_smooth_value_gradient(Z,fixed,beta)
        if mx<best:best=mx;bestX=Z.copy()
        grad=(grad-Z*np.sum(Z*grad,axis=1)[:,None])/norm[:,None]
        return value,grad.ravel()
    for trial in range(attempts):
        if time.monotonic()>=deadline:break
        y=initial*w if trial==0 else initial*w+rng.normal(scale=.05,size=initial.shape)
        try:
            for beta in betas:
                r=minimize(lambda v:evaluate(v,beta),y.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':maxiter,'ftol':1e-12,'gtol':1e-9});y=r.x.reshape(m,d)
            finished+=1
        except Cutoff:break
    if bestX is None:replacements=[]
    else:
        Q=np.rint(bestX/w*scale)
        if not np.all(np.isfinite(Q)) or np.max(np.abs(Q),initial=0)>=2**60:raise ValueError('rationalization out of range')
        replacements=Q.astype(np.int64).tolist()
    r={'parameters':p,'replacement_vectors':replacements,'best_max_correlation':None if bestX is None else best,'finished_attempts':finished,'forecast_cardinality':total,'status':'numerical proposal only'}
    return Artifact('search_result',{'backend':'cavity_repacking','result':r},{},Verification('unverified','cavity proposal; fixed-fixed and all other distances require exact promotion',{}))


def verify_search(a,context=None):
    context=context or RunContext(Path('.'));r=a.data.get('result',{});p=r.get('parameters');source,metric,removed,initial,d,total=_source(p,context);v=r.get('replacement_vectors')
    if not isinstance(v,list) or (v and (len(v)!=len(initial) or any(not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row) for row in v))):raise ValueError('invalid cavity result')
    return Verification('unverified','cavity result replayed, no distance claim',{'facts':{}})


def verify_exact(a,context=None):
    context=context or RunContext(Path('.'));data=a.data
    if a.kind!='implicit_spherical_code' or not isinstance(data,dict) or set(data)!={'backend','source','metric','remove_indices','replacement_vectors'} or data['backend']!='exact_cavity_repacking':raise ValueError('invalid exact cavity payload')
    source=data['source'];rm=data['remove_indices'];rep=data['replacement_vectors'];metric=data['metric']
    p={'source':source,'remove_indices':rm,'initial_replacements':rep,'metric':metric};_,_,removed,_,d,total=_source(p,context)
    if any(any(type(x) is not int for x in row) for row in rep):raise ValueError('exact replacement coordinates must be integer')
    rows=[v for i,v in enumerate(source) if i not in removed]+rep;mx=max(abs(x) for v in rows for x in v);context.limit('max_rational_bits',mx.bit_length())
    if not mx or d*max(metric)*mx*mx>=2**62:raise ValueError('exact integer dot products exceed safe range')
    V=np.asarray(rows,dtype=np.int64);weights=np.asarray(metric,dtype=np.int64);norms=[sum(w*x*x for w,x in zip(metric,row)) for row in rows]
    if any(n==0 for n in norms):raise ValueError('zero exact vector')
    distinct=sorted(set(norms));lookup={v:i for i,v in enumerate(distinct)};which=np.asarray([lookup[v] for v in norms],dtype=np.int32);context.limit('max_gram_entries',len(distinct)**2);T=np.array([[math.isqrt(a*b)//2 for b in distinct] for a in distinct],dtype=np.int64)
    batch=max(1,min(128,context.budget['max_gram_entries']//total));context.limit('max_gram_entries',batch*total)
    for at in range(0,total,batch):
        D=(V[at:at+batch]*weights)@V.T;threshold=T[which[at:at+batch,None],which[None,:]];bad=D>threshold
        for k in range(len(bad)):bad[k,:at+k+1]=False
        if np.any(bad):raise ValueError('exact cavity separation fails')
    rows_for_rank=rows[:max(64,4*d)];rank,_=exact_integer_rank(rows_for_rank,context)
    if rank<d:rank,_=exact_integer_rank(rows,context)
    facts={'dimension':rank,'rank_exact':rank,'ambient_dimension':d,'cardinality':total,'max_inner_product_upper_bound':'1/2'}
    return Verification('exact','weighted integer coordinates with exact pairwise thresholds',{'facts':facts,'data_sha256':digest(data)})


def certify_cavity(inputs,params,context):
    check_params(params,())
    if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='cavity_repacking':raise ValueError('expected cavity proposal')
    verify_search(inputs[0],context);r=inputs[0].data['result'];p=r['parameters'];v=r['replacement_vectors']
    if not v:raise ValueError('no replacement candidates')
    data={'backend':'exact_cavity_repacking','source':p['source'],'metric':p.get('metric',[1]*len(p['source'][0])),'remove_indices':p['remove_indices'],'replacement_vectors':v};a=Artifact('implicit_spherical_code',data,{},Verification());a.verification=verify_exact(a,context);a.facts=a.verification.details['facts'];return a

def adaptive_cavity(inputs,p,context):
    """Propose a point, release blockers, and expand the moving neighborhood."""
    if inputs:raise ValueError('adaptive cavity search takes no inputs')
    check_params(p,('source','metric','seed_candidate','expansion_rounds','add_per_round','attempts','max_iterations','seed','seconds','round_seconds','scale','betas'))
    source=p.get('source');candidate=p.get('seed_candidate')
    if not isinstance(source,list) or not source or not isinstance(source[0],list) or not isinstance(candidate,list):raise ValueError('source and seed_candidate required')
    d=len(source[0]);metric=p.get('metric',[1]*d)
    # Reuse source validation and allocation budgets before making arrays.
    check={'source':source,'metric':metric,'remove_indices':[],'initial_replacements':[candidate]}
    S,weights,_,seed,_,_=_source(check,context);w=np.sqrt(weights);Y=S*w;norm=np.linalg.norm(Y,axis=1);q=seed[0]*w;qn=np.linalg.norm(q)
    if np.any(norm<1e-14) or qn<1e-14:raise ValueError('zero source or candidate')
    Y/=norm[:,None];q/=qn;initial_scores=Y@q;removed=set(map(int,np.flatnonzero(initial_scores>.5+1e-12)))
    rounds=integer(p.get('expansion_rounds',2),'expansion_rounds');add=integer(p.get('add_per_round',32),'add_per_round');attempts=integer(p.get('attempts',1),'attempts');context.limit('max_trials',rounds*attempts)
    seconds=p.get('seconds',30);per=p.get('round_seconds',seconds)
    if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0 or type(per) not in (int,float) or not math.isfinite(per) or per<=0:raise ValueError('invalid adaptive time limits')
    context.limit('solver_seconds',seconds);deadline=time.monotonic()+seconds;initial=[source[i] for i in sorted(removed)]+[candidate];best=None;history=[];last=None
    for turn in range(rounds):
        remaining=deadline-time.monotonic()
        if remaining<=0:break
        opts={'source':source,'metric':metric,'remove_indices':sorted(removed),'initial_replacements':initial,'attempts':attempts,'max_iterations':p.get('max_iterations',100),'seconds':min(per,remaining),'seed':p.get('seed',0),'scale':p.get('scale',10000)}
        if 'betas' in p:opts['betas']=p['betas']
        a=search_cavity([],opts,context);r=a.data['result'];last=a;score=r['best_max_correlation'];history.append({'round':turn,'removed':len(removed),'moving':len(initial),'best_max_correlation':score})
        if score is not None and (best is None or score<best.data['result']['best_max_correlation']):best=a
        if turn+1>=rounds or not r['replacement_vectors']:break
        fixed=[i for i in range(len(source)) if i not in removed]
        if not fixed:break
        V=np.asarray(r['replacement_vectors'],float)*w;V/=np.linalg.norm(V,axis=1)[:,None];scores=np.max(V@Y[fixed].T,axis=0);take=np.argsort(-scores)[:min(add,len(fixed))];chosen=[fixed[int(i)] for i in take]
        removed.update(chosen);initial=r['replacement_vectors']+[source[i] for i in chosen]
    out=best or last
    if out is None:raise ValueError('adaptive search timed out before first round')
    out.data['result']['adaptive_history']=history
    return out

OPERATORS={'search.cavity_repacking':search_cavity,'search.adaptive_cavity':adaptive_cavity,'geometry.certify_cavity_repacking':certify_cavity}
INPUT_KINDS={'search.cavity_repacking':(),'search.adaptive_cavity':(),'geometry.certify_cavity_repacking':('search_result',)}
OUTPUT_KINDS={'search.cavity_repacking':('search_result',),'search.adaptive_cavity':('search_result',),'geometry.certify_cavity_repacking':('implicit_spherical_code',)}
