"""Search and exactly certify partially symmetric codes over Q(sqrt(3))."""
from pathlib import Path
import math
import numpy as np
from ..core import Artifact,RunContext,Verification,integer,check_params,digest
from .partial_group_assemblies import _params,_group,search_partial


def _positive(p,q):
    """Exact sign of p+q*sqrt(3), for Python integers."""
    if q==0:return p>0
    if q>0:return p>=0 or 3*q*q>p*p
    return p>0 and p*p>3*q*q


def _positive_array(p,q):
    return ((q==0)&(p>0)) | ((q>0)&((p>=0)|((p<0)&(3*q*q>p*p)))) | ((q<0)&(p>0)&(p*p>3*q*q))


def _source(p,context):
    check_params(p,('dimension','metric','orbit_specs','initial_a','initial_b','attempts','max_iterations','seed','seconds','scale','keep_initial','betas'))
    d=integer(p.get('dimension'),'dimension');A=p.get('initial_a');B=p.get('initial_b');spec=p.get('orbit_specs')
    if not isinstance(spec,list) or not spec or not isinstance(A,list) or not isinstance(B,list) or len(A)!=len(spec) or len(B)!=len(spec):raise ValueError('algebraic seed arrays must match orbit_specs')
    for R in (A,B):
        if any(not isinstance(row,list) or len(row)!=d or any(type(v) is not int for v in row) for row in R):raise ValueError('algebraic seeds must be integer matrices')
    keep=p.get('keep_initial',False)
    if type(keep) is not bool:raise ValueError('keep_initial must be boolean')
    base={k:p[k] for k in ('dimension','metric','orbit_specs','attempts','max_iterations','seed','seconds','scale','betas') if k in p}
    base['initial_seeds']=[(np.array(a,float)+math.sqrt(3)*np.array(b,float)).tolist() for a,b in zip(A,B)];base['keep_initial']=False
    d,spec,parts,N,metric=_params(base,context)
    for a,b,(_,den,H,basis,_) in zip(A,B,parts):
        for row in (a,b):
            if max(map(abs,row),default=0)*d*den>=2**62:raise ValueError('algebraic seed too large')
            v=np.asarray(row,dtype=np.int64)
            if not np.all(np.einsum('gij,j->gi',H,v)==den*v):raise ValueError('algebraic seed is not fixed by stabilizer')
        if not any(a) and not any(b):raise ValueError('zero algebraic seed')
    return base,A,B,parts,metric


def search_quadratic(inputs,p,context):
    if inputs:raise ValueError('quadratic group search takes no inputs')
    base,A,B,parts,metric=_source(p,context);r=search_partial(base,context);u=r['seed_vectors'];AA=[];BB=[]
    if u:
        if p.get('keep_initial',False) and r['best_max_correlation']>=.5-1e-8:AA=A;BB=B
        else:
            root=math.sqrt(3);w=np.asarray(metric,float)
            for a,b,v,(_,den,H,basis,reps) in zip(A,B,u,parts):
                real=np.asarray(a,float)+root*np.asarray(b,float);target=np.asarray(v,float)
                ratio=np.sqrt(np.sum(w*target*target)/np.sum(w*real*real));k=int(round(ratio));bb=np.asarray(b,dtype=np.int64)*k
                coeff=np.rint(np.linalg.lstsq(basis.astype(float),target-root*bb,rcond=None)[0])
                if not np.all(np.isfinite(coeff)) or float(np.max(np.abs(coeff),initial=0))*int(np.max(np.abs(basis),initial=0))*len(coeff)>=2**60:raise ValueError('algebraic rounding exceeds safe range')
                aa=coeff.astype(np.int64)@basis.T;AA.append(aa.tolist());BB.append(bb.tolist())
    result={'parameters':p,'seed_a':AA,'seed_b':BB,'best_max_correlation':r.get('best_max_correlation'),'status':'numerical proposal only'}
    return Artifact('search_result',{'backend':'quadratic_group_search','result':result},{},Verification('unverified','quadratic-field search; exact promotion required',{}))


def verify_search(a,context=None):
    context=context or RunContext(Path('.'));r=a.data.get('result',{});p=r.get('parameters');base,A,B,parts,metric=_source(p,context);sa=r.get('seed_a');sb=r.get('seed_b')
    if not isinstance(sa,list) or not isinstance(sb,list) or len(sa)!=len(sb) or (sa and (len(sa)!=len(parts) or any(not isinstance(v,list) or len(v)!=p['dimension'] or any(type(t) is not int for t in v) for v in sa+sb))):raise ValueError('invalid quadratic proposal')
    return Verification('unverified','quadratic-field proposal; exact separation not yet checked',{'facts':{}})


def verify_exact(a,context=None):
    context=context or RunContext(Path('.'));data=a.data
    if a.kind!='implicit_spherical_code' or not isinstance(data,dict) or set(data)!={'backend','groups','seed_a','seed_b','metric'} or data['backend']!='exact_quadratic_group_orbits':raise ValueError('invalid quadratic-field payload')
    groups,AA,BB,metric=data['groups'],data['seed_a'],data['seed_b'],data['metric']
    if not isinstance(AA,list) or not AA or not isinstance(BB,list) or not isinstance(groups,list) or len(groups)!=len(AA) or len(BB)!=len(AA) or not isinstance(AA[0],list) or not AA[0]:raise ValueError('missing quadratic orbit data')
    d=len(AA[0]);context.limit('max_dimension',d);context.limit('max_nodes',len(AA))
    if not isinstance(metric,list) or len(metric)!=d or any(type(w) is not int or w<=0 for w in metric):raise ValueError('invalid quadratic metric')
    if max(metric)*d>=2**62:raise ValueError('metric too large')
    W=np.asarray(metric,dtype=np.int64);diag=np.diag(W);pts=[];count=0
    for gens,aa,bb in zip(groups,AA,BB):
        if any(not isinstance(row,list) or len(row)!=d or any(type(v) is not int for v in row) for row in (aa,bb)):raise ValueError('invalid quadratic integer seed')
        M,den=_group(gens,d,context)
        if max(metric)*d*den*den>=2**62 or not np.all(np.einsum('gji,jk,gkl->gil',M,diag,M)==den*den*diag):raise ValueError('group does not preserve metric')
        mx=max(map(abs,aa+bb),default=0);context.limit('max_rational_bits',mx.bit_length())
        if not mx or mx*d*den>=2**62:raise ValueError('quadratic seed multiplication out of range')
        A=np.einsum('gij,j->gi',M,np.asarray(aa,dtype=np.int64));B=np.einsum('gij,j->gi',M,np.asarray(bb,dtype=np.int64));Z=np.unique(np.c_[A,B],axis=0);pts.append(Z);count+=len(Z);context.limit('max_materialized_points',count)
    context.limit('max_pairs',count*(count-1)//2);Z=np.concatenate(pts);A=Z[:,:d];B=Z[:,d:];ma=int(np.max(np.abs(A),initial=0));mb=int(np.max(np.abs(B),initial=0));wmax=max(metric)
    bound=d*wmax*(ma*ma+3*mb*mb+2*ma*mb)
    if 3*(32*bound*bound)**2>=2**62:
        # Python integers avoid overflow for larger exact numerators. This path
        # is deliberately budgeted because it does not use vectorized squares.
        context.limit('max_rank_work',count*(count-1)//2*d*4)
        rowsA=[[int(v) for v in row] for row in A];rowsB=[[int(v) for v in row] for row in B]
        norms=[]
        for av,bv in zip(rowsA,rowsB):
            npart=sum(w*(x*x+3*y*y) for w,x,y in zip(metric,av,bv));nrad=2*sum(w*x*y for w,x,y in zip(metric,av,bv))
            if not _positive(npart,nrad):raise ValueError('zero or negative quadratic norm')
            norms.append((npart,nrad))
        for i in range(count):
            for j in range(i+1,count):
                P=sum(w*(x*u+3*y*v) for w,x,y,u,v in zip(metric,rowsA[i],rowsB[i],rowsA[j],rowsB[j]));Q=sum(w*(x*v+y*u) for w,x,y,u,v in zip(metric,rowsA[i],rowsB[i],rowsA[j],rowsB[j]))
                if not _positive(P,Q):continue
                ni,ki=norms[i];nj,kj=norms[j];R=4*(P*P+3*Q*Q)-(ni*nj+3*ki*kj);S=8*P*Q-(ni*kj+ki*nj)
                if _positive(R,S):raise ValueError('exact quadratic-field separation fails')
    else:
        npart=np.sum((A*A+3*B*B)*W,axis=1);nrad=2*np.sum(A*B*W,axis=1)
        if not np.all(_positive_array(npart,nrad)):raise ValueError('zero or negative quadratic norm')
        batch=max(1,min(128,context.budget['max_gram_entries']//max(1,count)));context.limit('max_gram_entries',batch*count)
        for start in range(0,count,batch):
            a0=A[start:start+batch];b0=B[start:start+batch]
            P=(a0*W)@A.T+3*(b0*W)@B.T;Q=(a0*W)@B.T+(b0*W)@A.T
            ni=npart[start:start+batch,None];ki=nrad[start:start+batch,None];nj=npart[None,:];kj=nrad[None,:]
            R=4*(P*P+3*Q*Q)-(ni*nj+3*ki*kj);S=8*P*Q-(ni*kj+ki*nj)
            bad=_positive_array(P,Q)&_positive_array(R,S)
            for k in range(len(bad)):bad[k,:start+k+1]=False
            if np.any(bad):raise ValueError('exact quadratic-field separation fails')
    facts={'dimension':d,'ambient_dimension':d,'cardinality':count,'max_inner_product_upper_bound':'1/2'}
    return Verification('exact','materialized quadratic-field orbits and exact integer sign tests',{'facts':facts,'orbit_sizes':[len(x) for x in pts],'data_sha256':digest(data)})


def certify_quadratic(inputs,params,context):
    check_params(params,())
    if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='quadratic_group_search':raise ValueError('expected quadratic-field proposal')
    verify_search(inputs[0],context);r=inputs[0].data['result'];p=r['parameters']
    if not r['seed_a']:raise ValueError('no proposed seeds')
    data={'backend':'exact_quadratic_group_orbits','groups':[v.get('generators',[]) for v in p['orbit_specs']],'seed_a':r['seed_a'],'seed_b':r['seed_b'],'metric':p.get('metric',[1]*p['dimension'])};a=Artifact('implicit_spherical_code',data,{},Verification());a.verification=verify_exact(a,context);a.facts=a.verification.details['facts'];return a

OPERATORS={'search.quadratic_group_assemblies':search_quadratic,'groups.certify_quadratic_assemblies':certify_quadratic}
INPUT_KINDS={'search.quadratic_group_assemblies':(),'groups.certify_quadratic_assemblies':('search_result',)}
OUTPUT_KINDS={'search.quadratic_group_assemblies':('search_result',),'groups.certify_quadratic_assemblies':('implicit_spherical_code',)}
