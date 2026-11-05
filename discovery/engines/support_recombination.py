"""Recombine entire constant-weight support families under coordinate maps.

For 8-subsets with pairwise overlap at most 4, the incidence constraints on
5-subsets are equivalent to the full pairwise conflict graph. A numerical
solution is a proposal until checked independently. Auxiliary spherical layers
are not implied by an abstract support-code certificate.
"""
from itertools import combinations
from collections import defaultdict
from math import isfinite
from pathlib import Path
import numpy as np
from ..core import Artifact,RunContext,Verification,check_params,digest,integer


def _spec(p,ctx):
 check_params(p,('dimension','source','protected_supports','permutations','seconds','minimum_gain'))
 d=integer(p.get('dimension'),'dimension',8);ctx.limit('max_dimension',d);S=p.get('source');E=p.get('protected_supports',[]);P=p.get('permutations',[])
 if not isinstance(S,list) or not S or not isinstance(E,list) or not isinstance(P,list):raise ValueError('invalid source, protected supports, or permutations')
 if any(type(x) is not int or x<0 or x.bit_length()>d or x.bit_count()!=8 for x in S+E):raise ValueError('all supports must be 8-subsets')
 if len(set(S))!=len(S):raise ValueError('duplicate source supports')
 if any(not isinstance(p,list) or len(p)!=d or any(type(x) is not int for x in p) or sorted(p)!=list(range(d)) for p in P):raise ValueError('permutations must be bijections of coordinates')
 ctx.limit('max_codewords',len(S)*(1+len(P)));ctx.limit('max_rank_work',len(S)*(1+len(P))*56+len(S)*(1+len(P))*len(E));ctx.limit('max_nodes',max(1,len(P)))
 gain=integer(p.get('minimum_gain',1),'minimum_gain');seconds=p.get('seconds',5)
 if type(seconds) not in (int,float) or not isfinite(seconds) or seconds<=0:raise ValueError('invalid time limit')
 ctx.limit('solver_seconds',seconds)
 _check_code(S,d,E,ctx)
 return d,S,E,P,gain,seconds


def _check_code(S,d,E,ctx):
 if not S or len(set(S))!=len(S):raise ValueError('empty or duplicate support code')
 ctx.limit('max_codewords',len(S));ctx.limit('max_rank_work',len(S)*(56+len(E)))
 seen=set()
 for s in S:
  if type(s) is not int or s<0 or s.bit_length()>d or s.bit_count()!=8:raise ValueError('invalid support')
  if any((s&t).bit_count()>4 for t in E):raise ValueError('support conflicts with protected support')
  pos=[j for j in range(d) if s>>j&1]
  for sub in combinations(pos,5):
   mask=sum(1<<j for j in sub)
   if mask in seen:raise ValueError('overlap exceeds four')
   seen.add(mask)
 return True


def search(inputs,params,context):
 from scipy.optimize import milp,Bounds,LinearConstraint
 from scipy.sparse import coo_matrix,vstack
 if inputs:raise ValueError('support recombination takes no inputs')
 d,S,E,P,gain,seconds=_spec(params,context);pool=set(S)
 for perm in P:
  for s in S:
   q=sum(1<<perm[j] for j in range(d) if s>>j&1)
   if all((q&t).bit_count()<=4 for t in E):pool.add(q)
 V=sorted(pool);ctx=context;ctx.limit('max_codewords',len(V));by=defaultdict(list)
 for i,s in enumerate(V):
  pos=[j for j in range(d) if s>>j&1]
  for t in combinations(pos,5):by[sum(1<<j for j in t)].append(i)
 cl=[v for v in by.values() if len(v)>1];rows=np.repeat(np.arange(len(cl)),[len(v) for v in cl]);cols=np.array([i for v in cl for i in v],np.int64);M=coo_matrix((np.ones(len(cols)),(rows,cols)),shape=(len(cl),len(V))).tocsr();M=vstack([M,coo_matrix(-np.ones((1,len(V))))],format='csr');up=np.r_[np.ones(len(cl)),-float(len(S)+gain)]
 result=milp(-np.ones(len(V)),integrality=np.ones(len(V)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,up),options={'time_limit':float(seconds),'presolve':False})
 selected=[] if result.x is None else [V[i] for i in np.flatnonzero(result.x>.5)]
 if len(selected)<len(S)+gain:selected=[]
 data={'backend':'support_recombination','dimension':d,'source':S,'protected_supports':E,'selected':selected,'solver_status':str(result.message),'minimum_gain':gain,'candidate_count':len(V)}
 return Artifact('search_result',data,{},Verification('unverified','support family proposal; exact replay required'))


def certify(inputs,params,context):
 check_params(params,())
 if len(inputs)!=1 or inputs[0].kind!='search_result':raise ValueError('support proposal required')
 p=inputs[0].data
 if set(p)!={'backend','dimension','source','protected_supports','selected','solver_status','minimum_gain','candidate_count'} or p['backend']!='support_recombination':raise ValueError('invalid support proposal')
 d=integer(p['dimension'],'dimension',8);context.limit('max_dimension',d);S=p['source'];E=p['protected_supports'];Q=p['selected']
 if not isinstance(S,list) or not isinstance(E,list) or not isinstance(Q,list) or not isinstance(p['solver_status'],str):raise ValueError('invalid proposal fields')
 if any(type(x) is not int or x<0 or x.bit_length()>d or x.bit_count()!=8 for x in E):raise ValueError('invalid protected support')
 gain=integer(p['minimum_gain'],'minimum_gain');integer(p['candidate_count'],'candidate_count');_check_code(S,d,E,context)
 if Q and len(Q)<len(S)+gain:raise ValueError('proposal gain below requested target')
 chosen=Q if Q else S;_check_code(chosen,d,E,context)
 data={'dimension':d,'supports':chosen,'protected_supports':E};facts={'dimension':d,'cardinality':len(chosen),'weight':8,'max_intersection_upper_bound':4}
 return Artifact('support_code',data,facts,Verification('exact','5-subset packing and protected-support intersections',{'facts':facts,'data_sha256':digest(data)}))


def verify_code(artifact,context):
 if artifact.kind!='support_code' or set(artifact.data)!={'dimension','supports','protected_supports'}:raise ValueError('invalid support code')
 p=artifact.data;d=integer(p['dimension'],'dimension',8);context.limit('max_dimension',d);E=p['protected_supports'];S=p['supports']
 if not isinstance(E,list) or not isinstance(S,list) or any(type(x) is not int or x<0 or x.bit_length()>d or x.bit_count()!=8 for x in E):raise ValueError('invalid supports')
 _check_code(S,d,E,context);facts={'dimension':d,'cardinality':len(S),'weight':8,'max_intersection_upper_bound':4}
 from ..certification import compare_facts
 compare_facts(artifact.facts,facts)
 return Verification('exact','5-subset packing and protected-support intersections',{'facts':facts,'data_sha256':digest(p)})


def verify_search(a,context):
 code=certify([a],{},context)
 return Verification('unverified','support proposal replayed; spherical layers and optimality are not established',{'candidate_code_facts':code.facts})

OPERATORS={'search.support_recombination':search,'codes.certify_support_recombination':certify}
INPUT_KINDS={'search.support_recombination':(),'codes.certify_support_recombination':('search_result',)}
OUTPUT_KINDS={'search.support_recombination':('search_result',),'codes.certify_support_recombination':('support_code',)}
