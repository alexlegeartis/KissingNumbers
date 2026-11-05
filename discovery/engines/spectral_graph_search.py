"""Search rational and real-quadratic eigenspace representations of finite symmetric graphs.

The graph supplies a permutation representation. Exact spectral projectors
construct candidate spherical Gram matrices; numerical independent-set search
is only a proposal, and the returned Gram code is independently certified.
"""
from itertools import combinations
from math import isfinite
from time import monotonic
import numpy as np
from ..core import Artifact,check_params,integer
from ..exact_gram import verify_rational_gram,rational_json
from .quadratic_projector import verify_artifact as verify_quadratic,_sign
from fractions import Fraction


def _integer_projectors(A):
    """Use numerical eigenvalues only as proposals; validate all matrices exactly."""
    import sympy as sp
    n=len(A);mat=np.asarray(A,dtype=np.int64);vals=np.linalg.eigvalsh(mat.astype(float));rounded=np.rint(vals)
    if np.max(np.abs(vals-rounded),initial=0)>1e-7:return None
    roots=sorted(set(map(int,rounded)));I=np.eye(n,dtype=np.int64);deg=int(np.max(np.sum(np.abs(mat),axis=1),initial=0));bound=1
    for mu in roots:bound*=deg+abs(mu)
    if bound>=2**61:return None
    annihilator=I.copy()
    for mu in roots:annihilator=annihilator@(mat-mu*I)
    if np.any(annihilator):return None
    result=[]
    for lam in roots:
        num=I.copy();den=1
        for mu in roots:
            if mu!=lam:num=num@(mat-mu*I);den*=lam-mu
        if den<0:num=-num;den=-den
        vmax=int(np.max(np.abs(num),initial=0))
        if n*vmax*vmax>=2**61 or n*vmax*max(1,deg)>=2**61 or abs(den)*vmax>=2**61:return None
        if not np.array_equal(num@num,den*num) or not np.array_equal(mat@num,lam*num):return None
        tr=int(np.trace(num))
        if tr%den:return None
        mult=tr//den
        if mult>0:result.append((sp.Integer(lam),mult,sp.Matrix(num)/den))
    return result if result else None

def search(inputs,params,context):
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    import sympy as sp
    check_params(params,('adjacency','max_dimension','max_sectors','include_antipodes','seconds'))
    if inputs:raise ValueError('spectral graph search takes no inputs')
    A=params.get('adjacency')
    if not isinstance(A,list) or not A or not all(isinstance(row,list) and len(row)==len(A) for row in A):raise ValueError('adjacency must be nonempty and square')
    n=len(A);context.limit('max_materialized_points',2*n);context.limit('max_gram_entries',4*n*n);context.limit('max_rank_work',n*n*n*10)
    if any(type(x) is not int or x not in (0,1) for row in A for x in row):raise ValueError('adjacency must be binary integers')
    if any(A[i][i] or any(A[i][j]!=A[j][i] for j in range(i)) for i in range(n)):raise ValueError('adjacency must be symmetric with zero diagonal')
    maxdim=params.get('max_dimension',context.budget['max_dimension']);integer(maxdim,'max_dimension');context.limit('max_dimension',maxdim)
    maxsec=params.get('max_sectors',2);integer(maxsec,'max_sectors')
    antip=params.get('include_antipodes',True)
    if type(antip) is not bool:raise ValueError('include_antipodes must be boolean')
    seconds=params.get('seconds',5)
    if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
    context.limit('solver_seconds',seconds);start=monotonic();M=sp.Matrix(A);fast=_integer_projectors(A)
    ev={e:mult for e,mult,_ in fast} if fast is not None else M.eigenvals()
    radicals=set()
    for e in ev:
        for z in sp.preorder_traversal(e):
            if z.is_Pow and z.exp==sp.Rational(1,2) and z.base.is_Integer and z.base>0:radicals.add(int(z.base))
    if len(radicals)>1:raise ValueError('only one real quadratic field is supported')
    radicand=next(iter(radicals)) if radicals else None
    def coeff(x):
        x=sp.expand(sp.radsimp(x))
        if radicand is None:
            if not x.is_Rational:raise ValueError('only integer or single real quadratic graph spectra are supported')
            return (Fraction(int(x.p),int(x.q)),Fraction(0))
        root=sp.sqrt(radicand);b=sp.simplify(x.coeff(root));a=sp.simplify(x-b*root)
        if not a.is_Rational or not b.is_Rational:raise ValueError('only integer or single real quadratic graph spectra are supported')
        return (Fraction(int(a.p),int(a.q)),Fraction(int(b.p),int(b.q)))
    for e in ev:coeff(e)
    eig=sorted(((e,int(mult)) for e,mult in ev.items()),key=lambda z:sp.default_sort_key(z[0]));context.limit('max_trials',sum(comb(len(eig),r) for r in range(1,min(maxsec,len(eig))+1)))
    I=sp.eye(n);projectors=[]
    if fast is not None:
        by_root={e:P for e,mult,P in fast};projectors=[by_root[e] for e,mult in eig]
    else:
        for lam,mult in eig:
            P=I
            for mu,_ in eig:
                if mu!=lam:P=P*(M-mu*I)/(lam-mu)
            if radicand is not None:P=P.applyfunc(sp.simplify)
            if (P*P-P).applyfunc(sp.simplify)!=sp.zeros(n) or (M*P-lam*P).applyfunc(sp.simplify)!=sp.zeros(n) or sp.simplify(sp.trace(P)-mult)!=0:raise ValueError('spectral projector failed exact checks')
            projectors.append(P)
    best=None
    for size in range(1,min(maxsec,len(eig))+1):
        for inds in combinations(range(len(eig)),size):
            if sum(eig[i][1] for i in inds)>maxdim:continue
            P=sum((projectors[i] for i in inds),sp.zeros(n));classes={}
            for i in range(n):
                if (P[i,i]>0) if radicand is None else (_sign(coeff(P[i,i]),radicand)>0):classes.setdefault(P[i,i],[]).append(i)
            for diagonal,vertices in classes.items():
                points=[(i,1) for i in vertices]+([(i,-1) for i in vertices] if antip else [])
                N=len(points);conf=[]
                for i,(u,s) in enumerate(points):
                    for j,(v,t) in enumerate(points[:i]):
                        value=s*t*P[u,v]
                        if (2*value>diagonal) if radicand is None else (_sign(coeff(2*value-diagonal),radicand)>0):conf.append((j,i))
                chosen=[]
                # Greedy exact feasible incumbent even if the MILP times out.
                neighbors=[set() for _ in points]
                for i,j in conf:neighbors[i].add(j);neighbors[j].add(i)
                for i in sorted(range(N),key=lambda x:len(neighbors[x])):
                    if not any(j in neighbors[i] for j in chosen):chosen.append(i)
                remaining=float(seconds)-(monotonic()-start)
                if conf and remaining>0:
                    ii=np.array([i for i,j in conf],int);jj=np.array([j for i,j in conf],int);mat=coo_matrix((np.ones(2*len(ii)),(np.repeat(np.arange(len(ii)),2),np.c_[ii,jj].ravel())),shape=(len(ii),N)).tocsr()
                    r=milp(-np.ones(N),integrality=np.ones(N),bounds=Bounds(0,1),constraints=LinearConstraint(mat,-np.inf,1),options={'time_limit':remaining,'presolve':True})
                    if r.x is not None:
                        candidate=np.flatnonzero(r.x>.5).tolist()
                        if len(candidate)>len(chosen) and not any(i in candidate and j in candidate for i,j in conf):chosen=candidate
                if not chosen:continue
                if radicand is None:
                    gram=[[rational_json(P[points[i][0],points[j][0]]*points[i][1]*points[j][1]/diagonal) for j in chosen] for i in chosen]
                    data={'gram':gram,'norm_squared':1};report=verify_rational_gram(data,context);artifact=Artifact('gram_spherical_code',data,report.details['facts'],report)
                else:
                    projection=[[[rational_json(a),rational_json(b)] for a,b in (coeff(P[i,j]) for j in range(n))] for i in range(n)]
                    data={'radicand':radicand,'projection':projection,'selected':[list(points[i]) for i in chosen]};artifact=Artifact('quadratic_projector_code',data);report=verify_quadratic(artifact,context);artifact.facts=report.details['facts'];artifact.verification=report
                if best is None or len(chosen)>best.facts['cardinality']:best=artifact
    if best is None:raise ValueError('no permitted nonzero eigenspace')
    return best

# Local integer combination count avoids factorial intermediates.
from math import comb
OPERATORS={'search.spectral_graph':search}
INPUT_KINDS={'search.spectral_graph':()}
OUTPUT_KINDS={'search.spectral_graph':('gram_spherical_code','quadratic_projector_code')}
