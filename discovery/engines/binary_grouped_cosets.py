"""Synthesize large nonlinear binary codes from groups of subcode cosets."""
from math import isfinite,ceil
from time import monotonic
import numpy as np
from ..core import Artifact,Verification,check_params,integer
from .binary_coset_discovery import _basis,_dual,_distances
from .binary_span_cosets import build as build_span


def search(inputs,params,context):
 check_params(params,('length','parent_generator','codimension','subgroup_dimension','minimum_distance','subgroup_trials','seconds','seed','max_candidates','scramble_steps'))
 if inputs:raise ValueError('grouped coset search takes no inputs')
 n=integer(params.get('length'),'length');context.limit('max_dimension',n);G=params.get('parent_generator')
 if not isinstance(G,list) or not G:raise ValueError('parent generator required')
 _basis(G,n);k=len(G);c=integer(params.get('codimension'),'codimension');s=integer(params.get('subgroup_dimension'),'subgroup_dimension',0);d=integer(params.get('minimum_distance'),'minimum_distance');trials=integer(params.get('subgroup_trials',50),'subgroup_trials');seed=integer(params.get('seed',0),'seed',0);limit=integer(params.get('max_candidates',1000),'max_candidates');steps=integer(params.get('scramble_steps',5*k),'scramble_steps',0)
 if c>=k or s>c or n-k+c>21 or n>64 or d>n:raise ValueError('unsupported parameters')
 seconds=params.get('seconds',5)
 if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('invalid seconds')
 context.limit('solver_seconds',seconds);context.limit('max_trials',trials);states=1<<(n-k+c);context.limit('max_codewords',states);context.limit('max_rank_work',n*states+trials*states*(1<<s));context.limit('max_gram_entries',limit*limit)
 rng=np.random.default_rng(seed);H=G[:]
 for _ in range(steps):
  i,j=map(int,rng.choice(k,2,replace=False));H[i]^=H[j]
 base=H[:k-c];K=_dual(base,n);deadline=monotonic()+seconds;D=_distances(K,n);Q=np.flatnonzero(D>=d).astype(np.uint32);group_size=1<<s;needed=ceil(((1<<c)+1)/group_size)-1;stats=[];proposal=None
 for trial in range(trials):
  if monotonic()>=deadline:break
  group=[0];eligible=Q.copy();gens=[]
  for _ in range(s):
   if not len(eligible):break
   u=int(rng.choice(eligible));gens.append(u);group += [v^u for v in group];g=np.asarray(group,np.uint32);context.limit('max_gram_entries',len(eligible)*len(g));eligible=eligible[np.all(D[eligible[:,None]^g[None,:]]>=d,axis=1)]
  if len(gens)<s:stats.append({'trial':trial,'status':'incomplete_subgroup'});continue
  context.limit('max_gram_entries',len(eligible)*len(group));reps=np.unique(np.min(eligible[:,None]^np.asarray(group,np.uint32)[None,:],axis=1));m=len(reps)
  if m>limit:stats.append({'trial':trial,'candidates':m,'status':'candidate_budget'});continue
  good=np.ones((m,m),bool)
  for v in group:good &= D[reps[:,None]^reps[None,:]^v]>=d
  np.fill_diagonal(good,False);alive=good.sum(1)>=needed-1
  while alive.any():
   A=good[alive][:,alive];keep=A.sum(1)>=needed-1
   if keep.all():break
   ids=np.flatnonzero(alive);alive[ids[~keep]]=False
  ids=np.flatnonzero(alive);A=good[np.ix_(ids,ids)];adj=[sum(1<<int(j) for j in np.flatnonzero(row)) for row in A];timed=[False]
  def visit(path,choices):
   if monotonic()>=deadline:timed[0]=True;return None
   if len(path)>=needed:return path
   if len(path)+choices.bit_count()<needed:return None
   while choices:
    b=choices&-choices;i=b.bit_length()-1;choices-=b
    if len(path)+1+choices.bit_count()<needed:return None
    q=visit(path+[i],choices&adj[i])
    if q is not None:return q
    if timed[0]:return None
   return None
  result=visit([], (1<<len(ids))-1);stats.append({'trial':trial,'candidates':m,'core':len(ids),'status':'found' if result is not None else 'timeout' if timed[0] else 'exhausted'})
  if result is not None:
   labels=[int(a^b) for a in [0]+[int(v) for v in reps[ids[np.array(result)]]] for b in group];proposal={'generators':base,'parity_masks':K,'syndrome_labels':labels};break
 data={'backend':'binary_grouped_coset_search','length':n,'minimum_distance':d,'parent_size':1<<k,'proposal':proposal,'trials':stats}
 return Artifact('search_result',data,{},Verification('unverified','grouped coset discovery; exact span enumeration required'))


def certify(inputs,params,context):
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result' or inputs[0].data.get('backend')!='binary_grouped_coset_search':raise ValueError('grouped coset proposal required')
 data=inputs[0].data;p=data.get('proposal')
 if not isinstance(p,dict) or set(p)!={'generators','parity_masks','syndrome_labels'}:raise ValueError('no grouped coset proposal')
 n=integer(data.get('length'),'length');d=integer(data.get('minimum_distance'),'minimum_distance');parent=integer(data.get('parent_size'),'parent_size');G=p['generators'];K=p['parity_masks'];L=p['syndrome_labels']
 if any(not isinstance(v,list) for v in (G,K,L)) or len(G)+len(K)!=n or any(type(v) is not int or v<0 or v.bit_length()>n for v in G+K):raise ValueError('invalid quotient representation')
 _basis(G,n);_basis(K,n)
 if any((a&b).bit_count()%2 for a in G for b in K):raise ValueError('parity check is not orthogonal to kernel')
 if any(type(x) is not int or x<0 or x.bit_length()>len(K) for x in L):raise ValueError('invalid syndrome label')
 piv={}
 for j in range(n):
  v=sum(((K[i]>>j)&1)<<i for i in range(len(K)));word=1<<j
  while v:
   b=v.bit_length()-1
   if b in piv:v^=piv[b][0];word^=piv[b][1]
   else:piv[b]=(v,word);break
 if len(piv)!=len(K):raise ValueError('invalid parity rank')
 reps=[]
 for label in L:
  v=label;word=0
  while v:
   b=v.bit_length()-1;v^=piv[b][0];word^=piv[b][1]
  reps.append(word)
 code=build_span([] ,{'length':n,'generators':G,'representatives':reps},context)
 if code.facts['cardinality']<=parent or code.facts['min_distance']<d:raise ValueError('proposal does not establish requested improvement')
 return code


def verify_search(a,context):
 p=a.data
 if set(p)!={'backend','length','minimum_distance','parent_size','proposal','trials'} or not isinstance(p['trials'],list):raise ValueError('invalid grouped search receipt')
 if p['proposal'] is not None:certify([a],{},context)
 return Verification('unverified','grouped search receipt; no optimality claim')

OPERATORS={'search.binary_grouped_cosets':search,'codes.certify_binary_grouped_cosets':certify}
INPUT_KINDS={'search.binary_grouped_cosets':(),'codes.certify_binary_grouped_cosets':('search_result',)}
OUTPUT_KINDS={'search.binary_grouped_cosets':('search_result',),'codes.certify_binary_grouped_cosets':('binary_code',)}
