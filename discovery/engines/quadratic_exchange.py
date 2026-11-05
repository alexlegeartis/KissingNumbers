"""Search pointwise exchanges over Q(sqrt(3)); certify independently."""
from math import isfinite
import numpy as np
from ..core import Artifact,Verification,check_params
from .signed_layers import verify_artifact as verify_quadratic


def _positive(p,q):
 return ((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q))


def _points(p,context):
 A=p.get('source_a');B=p.get('source_b');C=p.get('candidate_a');D=p.get('candidate_b')
 if not all(isinstance(x,list) and x for x in (A,B,C,D)) or len(A)!=len(B) or len(C)!=len(D) or not isinstance(A[0],list) or not A[0]:raise ValueError('nonempty paired coefficient matrices required')
 n=len(A);m=len(C);d=len(A[0]);context.limit('max_dimension',d);context.limit('max_materialized_points',n+m);context.limit('max_codewords',m)
 if any(not isinstance(v,list) or len(v)!=d or any(type(x) is not int for x in v) for v in A+B+C+D):raise ValueError('integer coefficient matrices must have common width')
 mx=max(abs(x) for row in A+B+C+D for x in row);bound=4*d*mx*mx
 if 3*(8*bound*bound)**2>=2**62:raise ValueError('exact vectorized comparisons exceed safe integer range')
 A,B,C,D=(np.asarray(v,dtype=np.int64) for v in (A,B,C,D))
 for a,b in ((A,B),(C,D)):
  if np.any(np.sum(a*b,axis=1)!=0) or np.any(np.sum(a*a+3*b*b,axis=1)<=0):raise ValueError('positive rational squared norm required')
 return A,B,C,D


def _conflicts(a,b,c,d,norm0,norm1):
 P=a@c.T+3*b@d.T;Q=a@d.T+b@c.T
 R=4*(P*P+3*Q*Q)-norm0[:,None]*norm1[None,:];S=8*P*Q
 return _positive(P,Q)&_positive(R,S)


def search(inputs,params,context):
 from scipy.optimize import milp, Bounds, LinearConstraint
 from scipy.sparse import coo_matrix,vstack
 check_params(params,('source_a','source_b','candidate_a','candidate_b','seconds','minimum_gain','decompose'))
 if inputs:raise ValueError('quadratic exchange takes no inputs')
 A,B,C,D=_points(params,context);n=len(A);m=len(C);dim=A.shape[1]
 context.limit('max_pairs',m*(m-1)//2+m*n);context.limit('max_rank_work',(m*(m-1)//2+m*n)*dim)
 seconds=params.get('seconds',5);gain=params.get('minimum_gain',1);decompose=params.get('decompose',False)
 if type(decompose) is not bool or (decompose and gain!=1):raise ValueError('component mode requires gain 1 and boolean decompose')
 if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
 if type(gain) is not int or gain<1:raise ValueError('minimum_gain must be positive')
 context.limit('solver_seconds',seconds)
 base=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':params['source_a'],'radical_coordinates':params['source_b']});verify_quadratic(base,context)
 nx=np.sum(A*A+3*B*B,axis=1);nq=np.sum(C*C+3*D*D,axis=1)
 batch=max(1,min(96,context.budget['max_gram_entries']//max(n,m)));context.limit('max_gram_entries',batch*max(n,m));bi=[];bj=[];ii=[];jj=[]
 for st in range(0,m,batch):
  e=_conflicts(C[st:st+batch],D[st:st+batch],A,B,nq[st:st+batch],nx);a,b=np.where(e);bi.append(a+st);bj.append(b)
  e=_conflicts(C[st:st+batch],D[st:st+batch],C,D,nq[st:st+batch],nq)
  for k in range(len(e)):e[k,:st+k+1]=False
  a,b=np.where(e);ii.append(a+st);jj.append(b)
 a=np.concatenate(bi);b=np.concatenate(bj);i=np.concatenate(ii);j=np.concatenate(jj);ne=len(a);e=len(i)
 if decompose:
  from .exchange_components import search_components
  pick,removed,status=search_components(m,n,a,b,i,j,seconds)
 else:
  rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.c_[a,m+b].ravel(),np.c_[i,j].ravel()];v=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)]
  matrix=coo_matrix((v,(rows,cols)),shape=(ne+e,m+n)).tocsr();obj=np.r_[-np.ones(m),np.ones(n)];matrix=vstack([matrix,coo_matrix(obj.reshape(1,-1))],format='csr');up=np.r_[np.zeros(ne),np.ones(e),-float(gain)]
  result=milp(obj,integrality=np.ones(m+n),bounds=Bounds(0,1),constraints=LinearConstraint(matrix,-np.inf,up),options={'time_limit':float(seconds),'presolve':True})
  pick=[];removed=[]
  if result.x is not None:
   pick=np.flatnonzero(result.x[:m]>.5).tolist();removed=np.flatnonzero(result.x[m:]>.5).tolist()
   if len(pick)-len(removed)<gain:pick=[];removed=[]
  status=str(result.message)
 data={'backend':'quadratic_exchange','source_a':params['source_a'],'source_b':params['source_b'],'candidate_a':params['candidate_a'],'candidate_b':params['candidate_b'],'selected_indices':pick,'removed_indices':removed,'minimum_gain':gain,'solver_status':status}
 return Artifact('search_result',data,{},Verification('unverified','numerical exchange proposal; independent exact replay required'))


def certify(inputs,params,context):
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result':raise ValueError('quadratic exchange proposal required')
 v=inputs[0].data
 if set(v)!={'backend','source_a','source_b','candidate_a','candidate_b','selected_indices','removed_indices','minimum_gain','solver_status'} or v['backend']!='quadratic_exchange':raise ValueError('invalid proposal')
 A,B,C,D=_points(v,context);n=len(A);m=len(C);pick=v['selected_indices'];rm=v['removed_indices']
 if not isinstance(pick,list) or not isinstance(rm,list) or any(type(i) is not int or not 0<=i<m for i in pick) or any(type(i) is not int or not 0<=i<n for i in rm) or len(set(pick))!=len(pick) or len(set(rm))!=len(rm):raise ValueError('invalid selection')
 if type(v['minimum_gain']) is not int or v['minimum_gain']<1 or not isinstance(v['solver_status'],str):raise ValueError('invalid metadata')
 if (pick or rm) and len(pick)-len(rm)<v['minimum_gain']:raise ValueError('gain below target')
 keep=np.ones(n,dtype=bool);keep[rm]=False
 out=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':np.r_[A[keep],C[pick]].tolist(),'radical_coordinates':np.r_[B[keep],D[pick]].tolist()})
 out.verification=verify_quadratic(out,context);out.facts=out.verification.details['facts'];return out


def verify_search(artifact,context):
 result=certify([artifact],{},context)
 return Verification('unverified','proposal replayed exactly; no solver optimality assertion',{'candidate_code_facts':result.facts})

OPERATORS={'search.quadratic_exchange':search,'geometry.certify_quadratic_exchange':certify}
INPUT_KINDS={'search.quadratic_exchange':(),'geometry.certify_quadratic_exchange':('search_result',)}
OUTPUT_KINDS={'search.quadratic_exchange':('search_result',),'geometry.certify_quadratic_exchange':('quadratic_spherical_code',)}
