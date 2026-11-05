"""Pointwise exchange search on exactly normalized integer directions.

Source and candidate vectors may have unequal norms. All proposed changes
are replayed by the exact weighted cavity verifier before an artifact is
returned; a solver timeout is not an upper-bound proof.
"""
from math import isfinite
import numpy as np
from ..core import Artifact,Verification,check_params,integer
from .cavity_repacking import verify_exact


def search(inputs,params,context):
 from scipy.optimize import milp,Bounds,LinearConstraint
 from scipy.sparse import coo_matrix,vstack
 check_params(params,('source','candidates','metric','seconds','minimum_gain'))
 if inputs:raise ValueError('projective exchange takes no inputs')
 X=params.get('source');Q=params.get('candidates');metric=params.get('metric')
 if not isinstance(X,list) or not X or not isinstance(X[0],list) or not X[0]:raise ValueError('nonempty source required')
 d=len(X[0]);context.limit('max_dimension',d)
 if not isinstance(Q,list) or not Q:raise ValueError('nonempty candidate pool required')
 if any(not isinstance(v,list) or len(v)!=d or any(type(x) is not int for x in v) for v in X+Q):raise ValueError('source and candidates require integer rows of common width')
 if metric is None:metric=[1]*d
 if not isinstance(metric,list) or len(metric)!=d or any(type(x) is not int or x<=0 for x in metric):raise ValueError('positive integer diagonal metric required')
 n=len(X);m=len(Q);context.limit('max_materialized_points',n+m);context.limit('max_codewords',m);context.limit('max_pairs',m*(m-1)//2+m*n);context.limit('max_rank_work',(m*(m-1)//2+m*n)*d)
 seconds=params.get('seconds',5);gain=params.get('minimum_gain',1)
 if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
 if type(gain) is not int or gain<1:raise ValueError('minimum_gain must be positive')
 context.limit('solver_seconds',seconds)
 # Independently verify the supplied source before any optimization.
 original={'backend':'exact_cavity_repacking','source':X,'metric':metric,'remove_indices':[0],'replacement_vectors':[X[0]]};baseline=Artifact('implicit_spherical_code',original);report=verify_exact(baseline,context);baseline.facts=report.details['facts'];baseline.verification=report
 mx=max(abs(x) for row in X+Q for x in row);bound=d*max(metric)*mx*mx
 if 4*bound*bound>=2**62:raise ValueError('exact squared comparisons exceed safe vectorized range')
 V=np.asarray(X,np.int64);U=np.asarray(Q,np.int64);W=np.asarray(metric,np.int64)
 NX=np.sum(V*V*W,axis=1);NQ=np.sum(U*U*W,axis=1)
 if np.any(NX<=0) or np.any(NQ<=0):raise ValueError('zero vector')
 bi=[];bj=[];ii=[];jj=[];batch=max(1,min(128,context.budget['max_gram_entries']//max(n,m)));context.limit('max_gram_entries',batch*max(n,m))
 for st in range(0,m,batch):
  D=(U[st:st+batch]*W)@V.T;bad=(D>0)&(4*D*D>NQ[st:st+batch,None]*NX[None,:]);a,b=np.where(bad);bi.append(a+st);bj.append(b)
  D=(U[st:st+batch]*W)@U.T;bad=(D>0)&(4*D*D>NQ[st:st+batch,None]*NQ[None,:])
  for k in range(len(bad)):bad[k,:st+k+1]=False
  a,b=np.where(bad);ii.append(a+st);jj.append(b)
 a=np.concatenate(bi);b=np.concatenate(bj);i=np.concatenate(ii);j=np.concatenate(jj);ne=len(a);e=len(i)
 rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.c_[a,m+b].ravel(),np.c_[i,j].ravel()];values=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];A=coo_matrix((values,(rows,cols)),shape=(ne+e,m+n)).tocsr();obj=np.r_[-np.ones(m),np.ones(n)];A=vstack([A,coo_matrix(obj.reshape(1,-1))],format='csr');up=np.r_[np.zeros(ne),np.ones(e),-float(gain)]
 r=milp(obj,integrality=np.ones(m+n),bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,up),options={'time_limit':float(seconds),'presolve':True})
 pick=[];removed=[]
 if r.x is not None:
  pick=np.flatnonzero(r.x[:m]>.5).tolist();removed=np.flatnonzero(r.x[m:]>.5).tolist()
  if len(pick)-len(removed)<gain:pick=[];removed=[]
 data={'backend':'projective_exchange','source':X,'candidates':Q,'metric':metric,'selected_indices':pick,'removed_indices':removed,'solver_status':str(r.message),'minimum_gain':gain}
 return Artifact('search_result',data,{},Verification('unverified','numerical exchange proposal; exact replay required'))


def certify(inputs,params,context):
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='projective_exchange':raise ValueError('projective exchange proposal required')
 data=inputs[0].data
 if set(data)!={'backend','source','candidates','metric','selected_indices','removed_indices','solver_status','minimum_gain'}:raise ValueError('invalid projective proposal')
 X=data['source'];Q=data['candidates'];pick=data['selected_indices'];rm=data['removed_indices']
 if not isinstance(X,list) or not X or not isinstance(Q,list) or not isinstance(pick,list) or any(type(i) is not int or i<0 or i>=len(Q) for i in pick) or len(set(pick))!=len(pick):raise ValueError('invalid candidate selection')
 if not isinstance(rm,list) or any(type(i) is not int or i<0 or i>=len(X) for i in rm) or len(set(rm))!=len(rm):raise ValueError('invalid source removals')
 if type(data['minimum_gain']) is not int or data['minimum_gain']<1 or not isinstance(data['solver_status'],str):raise ValueError('invalid search metadata')
 if (pick or rm) and len(pick)-len(rm)<data['minimum_gain']:raise ValueError('proposal does not meet its target gain')
 if not pick and not rm:rm=[0];replacement=[X[0]]
 else:replacement=[Q[i] for i in pick]
 if not replacement:raise ValueError('empty replacement with nonempty removals')
 output={'backend':'exact_cavity_repacking','source':X,'metric':data['metric'],'remove_indices':rm,'replacement_vectors':replacement};result=Artifact('implicit_spherical_code',output,{},Verification());result.verification=verify_exact(result,context);result.facts=result.verification.details['facts'];return result

def verify_search(artifact,context):
 code=certify([artifact],{},context)
 return Verification('unverified','proposal replayed exactly; no solver optimality assertion',{'candidate_code_facts':code.facts})

OPERATORS={'search.projective_exchange':search,'geometry.certify_projective_exchange':certify}
INPUT_KINDS={'search.projective_exchange':(),'geometry.certify_projective_exchange':('search_result',)}
OUTPUT_KINDS={'search.projective_exchange':('search_result',),'geometry.certify_projective_exchange':('implicit_spherical_code',)}
