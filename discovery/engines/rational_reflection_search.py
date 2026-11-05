"""Generate rational reflection candidates and search exact projective exchanges.

A reflection of x in the hyperplane orthogonal to e has numerator
(e.g.e)x - 2(x.g.e)e.  It is an exact rational orthogonal transformation;
integer gcd reduction retains the same oriented ray.  Candidate selection and
MILP search are finite and do not imply global optimality.
"""
from math import gcd
import numpy as np
from ..core import Artifact,Verification,check_params,integer
from .projective_exchange import search as exchange_search


def search(inputs,params,context):
 check_params(params,('source','axes','metric','candidate_limit','seconds','minimum_gain','decompose'))
 if inputs:raise ValueError('rational reflection search takes no inputs')
 X=params.get('source');axes=params.get('axes')
 if not isinstance(X,list) or not X or not isinstance(X[0],list) or not X[0]:raise ValueError('nonempty source required')
 d=len(X[0]);context.limit('max_dimension',d)
 if not isinstance(axes,list) or not axes:raise ValueError('nonempty axes required')
 if any(not isinstance(v,list) or len(v)!=d or any(type(a) is not int for a in v) for v in X+axes):raise ValueError('integer vectors of common width required')
 g=params.get('metric',[1]*d)
 if not isinstance(g,list) or len(g)!=d or any(type(a) is not int or a<=0 for a in g):raise ValueError('positive integer diagonal metric required')
 cap=integer(params.get('candidate_limit',200),'candidate_limit');context.limit('max_trials',len(axes));context.limit('max_materialized_points',len(X)+len(axes)*len(X));context.limit('max_rank_work',len(axes)*len(X)*d)
 def primitive(v):
  k=0
  for x in v:k=gcd(k,abs(x))
  if not k:raise ValueError('zero vector')
  return tuple(x//k for x in v)
 source=[primitive(v) for v in X]
 if len(set(source))!=len(source):raise ValueError('duplicate source directions')
 found=set();old=set(source)
 for e in axes:
  e=primitive(e);en=sum(g[j]*e[j]*e[j] for j in range(d))
  for x in source:
   dot=sum(g[j]*x[j]*e[j] for j in range(d));v=primitive([en*x[j]-2*dot*e[j] for j in range(d)])
   if v not in old:found.add(v)
 if not found:raise ValueError('reflection axes produced no new directions')
 candidates=sorted(found);context.limit('max_codewords',len(candidates));context.limit('max_pairs',len(candidates)*len(source));context.limit('max_rank_work',len(candidates)*len(source)*d)
 mx=max(abs(v) for row in source+candidates for v in row);bound=d*max(g)*mx*mx
 if 4*bound*bound>=2**62:raise ValueError('reflected coordinates exceed safe exact comparison range')
 Xn=np.asarray(source,np.int64);Q=np.asarray(candidates,np.int64);W=np.asarray(g,np.int64);NX=((Xn*Xn)*W).sum(1);NQ=((Q*Q)*W).sum(1);cost=np.empty(len(Q),np.int64)
 batch=max(1,min(64,context.budget['max_gram_entries']//max(len(Xn),1)));context.limit('max_gram_entries',batch*len(Xn))
 for st in range(0,len(Q),batch):
  D=(Q[st:st+batch]*W)@Xn.T;cost[st:st+batch]=((D>0)&(4*D*D>NQ[st:st+batch,None]*NX[None,:])).sum(1)
 chosen=np.argsort(cost,kind='stable')[:cap];selected=[candidates[int(i)] for i in chosen]
 p={'source':[list(v) for v in source],'candidates':[list(v) for v in selected],'metric':g,'seconds':params.get('seconds',5),'minimum_gain':params.get('minimum_gain',1),'decompose':params.get('decompose',False)}
 result=exchange_search([],p,context);result.data['solver_status']='rational reflection pool; '+result.data['solver_status'];return result

def quadratic_search(inputs,params,context):
 from .quadratic_exchange import search as quadratic_exchange_search,_conflicts
 check_params(params,('source_a','source_b','axes_a','axes_b','metric','candidate_limit','seconds','minimum_gain','decompose'))
 if inputs:raise ValueError('quadratic reflection search takes no inputs')
 A=params.get('source_a');B=params.get('source_b');EA=params.get('axes_a');EB=params.get('axes_b')
 if not all(isinstance(x,list) and x for x in (A,B,EA,EB)) or len(A)!=len(B) or len(EA)!=len(EB) or not isinstance(A[0],list) or not A[0]:raise ValueError('paired source and axis matrices required')
 d=len(A[0]);context.limit('max_dimension',d)
 if any(not isinstance(v,list) or len(v)!=d or any(type(a) is not int for a in v) for v in A+B+EA+EB):raise ValueError('integer coefficient rows of common width required')
 metric=params.get('metric',[1]*d)
 if not isinstance(metric,list) or len(metric)!=d or any(type(x) is not int or x<=0 for x in metric):raise ValueError('positive integer diagonal metric required')
 expanded=sum(metric);context.limit('max_dimension',expanded)
 if expanded!=d:
  def expand(rows):return [[x for x,k in zip(row,metric) for _ in range(k)] for row in rows]
  A,B,EA,EB=map(expand,(A,B,EA,EB));d=expanded
 cap=integer(params.get('candidate_limit',200),'candidate_limit');context.limit('max_trials',len(EA));context.limit('max_materialized_points',len(A)+len(EA)*len(A));context.limit('max_rank_work',len(EA)*len(A)*d)
 def primitive(a,b):
  k=0
  for x in a+b:k=gcd(k,abs(x))
  if not k:raise ValueError('zero vector')
  return tuple(x//k for x in a),tuple(x//k for x in b)
 source=[primitive(a,b) for a,b in zip(A,B)];axes=[primitive(a,b) for a,b in zip(EA,EB)]
 if len(set(source))!=len(source):raise ValueError('duplicate source directions')
 if any(sum(a[i]*b[i] for i in range(d))!=0 for a,b in source+axes):raise ValueError('rational squared norms required')
 old=set(source);found=set()
 for e,f in axes:
  en=sum(e[j]*e[j]+3*f[j]*f[j] for j in range(d))
  for a,b in source:
   p=sum(a[j]*e[j]+3*b[j]*f[j] for j in range(d));q=sum(a[j]*f[j]+b[j]*e[j] for j in range(d))
   ra=[en*a[j]-2*(p*e[j]+3*q*f[j]) for j in range(d)];rb=[en*b[j]-2*(p*f[j]+q*e[j]) for j in range(d)];v=primitive(ra,rb)
   if v not in old:found.add(v)
 if not found:raise ValueError('reflection axes produced no new directions')
 candidates=sorted(found);context.limit('max_codewords',len(candidates));context.limit('max_pairs',len(candidates)*len(source));context.limit('max_rank_work',len(candidates)*len(source)*d)
 norms=[sum(a[j]*a[j]+3*b[j]*b[j] for j in range(d)) for a,b in source+candidates]
 if 243*max(norms)**4>=2**62:raise ValueError('quadratic reflections exceed safe exact comparison range')
 A0=np.asarray([v[0] for v in source],np.int64);B0=np.asarray([v[1] for v in source],np.int64);C=np.asarray([v[0] for v in candidates],np.int64);D=np.asarray([v[1] for v in candidates],np.int64);NX=(A0*A0+3*B0*B0).sum(1);NQ=(C*C+3*D*D).sum(1);cost=np.empty(len(C),np.int64)
 batch=max(1,min(64,context.budget['max_gram_entries']//max(len(A0),1)));context.limit('max_gram_entries',batch*len(A0))
 for st in range(0,len(C),batch):cost[st:st+batch]=_conflicts(C[st:st+batch],D[st:st+batch],A0,B0,NQ[st:st+batch],NX).sum(1)
 ids=np.argsort(cost,kind='stable')[:cap];p={'source_a':A0.tolist(),'source_b':B0.tolist(),'candidate_a':C[ids].tolist(),'candidate_b':D[ids].tolist(),'seconds':params.get('seconds',5),'minimum_gain':params.get('minimum_gain',1),'decompose':params.get('decompose',False)}
 result=quadratic_exchange_search([],p,context);result.data['solver_status']='quadratic reflection pool; '+result.data['solver_status']
 if expanded!=len(metric):return Artifact('search_result',{'backend':'weighted_quadratic_reflection','proposal':result.data,'metric':metric},{},Verification('unverified','weighted quadratic reflection proposal'))
 return result

def verify_weighted(artifact,context):
 from .signed_layers import verify_artifact as verify_quadratic,_full_rank_witness
 from ..certification import compare_facts
 if artifact.kind!='weighted_quadratic_spherical_code':raise ValueError('weighted quadratic code required')
 data=artifact.data
 if set(data)!={'radicand','rational_coordinates','radical_coordinates','metric'} or data['radicand']!=3:raise ValueError('invalid weighted quadratic code')
 A=data['rational_coordinates'];B=data['radical_coordinates'];metric=data['metric']
 if not isinstance(A,list) or not A or not isinstance(A[0],list) or not A[0] or not isinstance(B,list) or len(A)!=len(B):raise ValueError('invalid coordinate matrices')
 d=len(A[0]);context.limit('max_dimension',d)
 if not isinstance(metric,list) or len(metric)!=d or any(type(x) is not int or x<=0 for x in metric):raise ValueError('invalid metric')
 if any(not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row) for row in A+B):raise ValueError('invalid coordinate rows')
 def expand(rows):return [[x for x,k in zip(row,metric) for _ in range(k)] for row in rows]
 tmp=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':expand(A),'radical_coordinates':expand(B)})
 report=verify_quadratic(tmp,context);rank=_full_rank_witness(A,B,context);facts={'dimension':d,'ambient_dimension':d,'cardinality':len(A),'max_cosine_upper_bound':[1,2],'rank_lower_bound':rank,'dimension_scope':'weighted ambient; modular minor lower bound'}
 if rank==d:facts.update(rank_exact=d,dimension_scope='exact; nonzero algebraic minor proved modulo13')
 norms=[sum(metric[j]*(a[j]*a[j]+3*b[j]*b[j]) for j in range(d)) for a,b in zip(A,B)]
 if len(set(norms))==1:facts['norm_squared']=norms[0]
 compare_facts(artifact.facts,facts)
 return Verification('exact','weighted quadratic code via exact isometric expansion',{'facts':facts,'pairs_checked':len(A)*(len(A)-1)//2,'expanded_dimension':sum(metric)})


def certify_weighted(inputs,params,context):
 from .quadratic_exchange import certify as certify_quadratic
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result':raise ValueError('weighted quadratic reflection proposal required')
 data=inputs[0].data
 if set(data)!={'backend','proposal','metric'} or data['backend']!='weighted_quadratic_reflection':raise ValueError('invalid proposal')
 metric=data['metric']
 if not isinstance(metric,list) or not metric or any(type(x) is not int or x<=0 for x in metric):raise ValueError('invalid metric')
 expanded=certify_quadratic([Artifact('search_result',data['proposal'])],{},context)
 def compress(rows):
  out=[]
  for row in rows:
   if len(row)!=sum(metric):raise ValueError('expanded dimension mismatch')
   pos=0;v=[]
   for k in metric:
    block=row[pos:pos+k]
    if any(x!=block[0] for x in block):raise ValueError('expanded coordinate blocks changed')
    v.append(block[0]);pos+=k
   out.append(v)
  return out
 result=Artifact('weighted_quadratic_spherical_code',{'radicand':3,'rational_coordinates':compress(expanded.data['rational_coordinates']),'radical_coordinates':compress(expanded.data['radical_coordinates']),'metric':metric})
 result.verification=verify_weighted(result,context);result.facts=result.verification.details['facts'];return result


def verify_search(artifact,context):
 result=certify_weighted([artifact],{},context)
 return Verification('unverified','weighted proposal replayed exactly; no optimality assertion',{'candidate_code_facts':result.facts})


OPERATORS={'search.rational_reflection_exchange':search,'search.quadratic_reflection_exchange':quadratic_search,'geometry.certify_weighted_quadratic_reflection':certify_weighted}
INPUT_KINDS={'search.rational_reflection_exchange':(),'search.quadratic_reflection_exchange':(),'geometry.certify_weighted_quadratic_reflection':('search_result',)}
OUTPUT_KINDS={'search.rational_reflection_exchange':('search_result',),'search.quadratic_reflection_exchange':('search_result',),'geometry.certify_weighted_quadratic_reflection':('weighted_quadratic_spherical_code',)}
