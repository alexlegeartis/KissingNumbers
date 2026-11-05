"""Discover larger binary codes as unions of subcode cosets.

This searches the entire syndrome quotient for each sampled subcode. A clique
in the exact syndrome-distance graph proposes a union; independent BFS and
coordinate DP certify any resulting code. It does not assume binary linearity
of the union or that a search failure proves nonexistence.
"""
from math import isfinite
from time import monotonic
import numpy as np
from ..core import Artifact,Verification,check_params,integer
from . import cosets


def _basis(rows,n):
 piv={}
 for x in rows:
  if type(x) is not int or x<0 or x.bit_length()>n:raise ValueError('invalid binary generator')
  v=x
  while v:
   j=v.bit_length()-1
   if j in piv:v^=piv[j]
   else:piv[j]=v;break
 if len(piv)!=len(rows):raise ValueError('binary generators must be independent')
 return piv


def _dual(rows,n):
 piv=_basis(rows,n);free=[j for j in range(n) if j not in piv];K=[]
 for j in free:
  v=1<<j
  for i in sorted(piv):
   if (v&piv[i]).bit_count()%2:v|=1<<i
  K.append(v)
 assert all((v&w).bit_count()%2==0 for v in K for w in rows)
 return K


def _distances(K,n):
 r=len(K);states=1<<r;cols=np.array([sum(((K[i]>>j)&1)<<i for i in range(r)) for j in range(n)],np.uint32);D=np.full(states,255,np.uint8);D[0]=0;front=np.array([0],np.uint32)
 for depth in range(1,n+1):
  q=np.unique((front[:,None]^cols[None,:]).ravel());front=q[D[q]==255];D[front]=depth
  if not len(front):break
 if np.any(D==255):raise ValueError('syndrome graph not connected')
 return D


def _clique(Q,D,target,deadline):
 n=len(Q);G=D[Q[:,None]^Q[None,:]]>=target;np.fill_diagonal(G,False)
 return G


def search(inputs,params,context):
 check_params(params,('length','parent_generator','codimension','minimum_distance','trials','seconds','seed','max_candidates','scramble_steps'))
 if inputs:raise ValueError('binary coset discovery takes no inputs')
 n=integer(params.get('length'),'length');context.limit('max_dimension',n);G=params.get('parent_generator')
 if not isinstance(G,list) or not G:raise ValueError('nonempty binary generator required')
 _basis(G,n);k=len(G);c=integer(params.get('codimension',3),'codimension');d=integer(params.get('minimum_distance'),'minimum_distance');trials=integer(params.get('trials',1),'trials');seed=integer(params.get('seed',0),'seed',0);limit=integer(params.get('max_candidates',1000),'max_candidates');steps=integer(params.get('scramble_steps',5*k),'scramble_steps',0)
 if c>=k or n-k+c>19 or d>n or n>64:raise ValueError('unsupported codimension, quotient rank, distance, or length')
 seconds=params.get('seconds',5)
 if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
 context.limit('solver_seconds',seconds);context.limit('max_trials',trials);states=1<<(n-k+c);context.limit('max_codewords',states);context.limit('max_rank_work',trials*n*states);context.limit('max_gram_entries',limit*limit)
 rng=np.random.default_rng(seed);deadline=monotonic()+seconds;stats=[];found=None;target=1<<c
 for trial in range(trials):
  if monotonic()>=deadline:break
  H=G[:]
  for _ in range(steps):
   i,j=map(int,rng.choice(k,2,replace=False));H[i]^=H[j]
  K=_dual(H[:k-c],n);D=_distances(K,n);Q=np.flatnonzero(D>=d).astype(np.uint32);num=len(Q)
  if num>limit:stats.append({'trial':trial,'candidate_count':num,'status':'candidate_budget'});continue
  good=_clique(Q,D,d,deadline);alive=good.sum(1)>=target-1
  while alive.any():
   B=good[alive][:,alive];keep=B.sum(1)>=target-1
   if keep.all():break
   ids=np.flatnonzero(alive);alive[ids[~keep]]=False
  ids=np.flatnonzero(alive);B=good[np.ix_(ids,ids)];adj=[sum(1<<int(j) for j in np.flatnonzero(row)) for row in B];timed=[False]
  def visit(path,choices):
   if monotonic()>=deadline:timed[0]=True;return None
   if len(path)>=target:return path
   if len(path)+choices.bit_count()<target:return None
   while choices:
    b=choices&-choices;i=b.bit_length()-1;choices-=b
    if len(path)+1+choices.bit_count()<target:return None
    q=visit(path+[i],choices&adj[i])
    if q is not None:return q
    if timed[0]:return None
   return None
  clique=visit([], (1<<len(ids))-1);stats.append({'trial':trial,'candidate_count':num,'core_count':len(ids),'status':'found' if clique else 'timeout' if timed[0] else 'exhausted'})
  if clique is not None:
   chosen=[0]+[int(v) for v in Q[ids[np.array(clique)]]];found={'parity_check':[[int((v>>j)&1) for j in range(n)] for v in K],'syndrome_labels':chosen};break
 data={'backend':'binary_subcode_coset_discovery','length':n,'minimum_distance':d,'parent_size':1<<k,'proposal':found,'trials':stats}
 return Artifact('search_result',data,{},Verification('unverified','binary subcode coset discovery; exact code certification required'))


def certify(inputs,params,context):
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='binary_subcode_coset_discovery':raise ValueError('binary coset proposal required')
 data=inputs[0].data;proposal=data.get('proposal')
 if not isinstance(proposal,dict) or set(proposal)!={'parity_check','syndrome_labels'}:raise ValueError('no coset proposal')
 n=integer(data.get('length'),'length');d=integer(data.get('minimum_distance'),'minimum_distance');parent=integer(data.get('parent_size'),'parent_size');oracle=cosets.build([] ,{'parity_check':proposal['parity_check']},context);code=cosets.union([oracle],{'syndrome_labels':proposal['syndrome_labels']},context)
 if code.facts['length']!=n or code.facts['cardinality']<=parent or code.facts['min_distance']<d:raise ValueError('proposal does not prove the requested improvement')
 return code


def verify_search(a,context):
 data=a.data
 if set(data)!={'backend','length','minimum_distance','parent_size','proposal','trials'} or not isinstance(data['trials'],list):raise ValueError('invalid discovery result')
 if data['proposal'] is not None:certify([a],{},context)
 return Verification('unverified','search receipt; no claim of global optimality')

OPERATORS={'search.binary_subcode_cosets':search,'codes.certify_binary_subcode_cosets':certify}
INPUT_KINDS={'search.binary_subcode_cosets':(),'codes.certify_binary_subcode_cosets':('search_result',)}
OUTPUT_KINDS={'search.binary_subcode_cosets':('search_result',),'codes.certify_binary_subcode_cosets':('binary_code',)}
