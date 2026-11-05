"""Construct nonlinear binary Kerdock codes via trace over GR(4,m).

The Z4-linear generator is synthesized from Teichmuller representatives and
Galois-ring trace. The Gray image is not assumed binary-linear. Exact Lee
weight enumeration certifies cardinality and the complete binary distance.
"""
from itertools import product
from pathlib import Path
import numpy as np
from ..core import Artifact,RunContext,Verification,check_params,digest
from ..certification import verify_integer_spherical
from .finite_field_search import _field


def _construct(m,context,extra_evaluations=None):
 if type(m) is not int or m not in (3,5):raise ValueError('supported odd extension degrees are 3 and 5')
 base={3:11,5:37}[m];q=_field(m,base);count=4**(m+1);extra=extra_evaluations if extra_evaluations is not None else []
 if not isinstance(extra,list) or any(not isinstance(v,list) or len(v)!=m or any(type(x) is not int or not 0<=x<4 for x in v) for v in extra) or len({tuple(v) for v in extra})!=len(extra):raise ValueError('extra evaluations must be distinct ring elements')
 context.limit('max_codewords',count);context.limit('max_dimension',2*(q+len(extra)));context.limit('max_rank_work',count*(q+len(extra))*(m+1))
 zero=(0,)*m;one=(1,)+(0,)*(m-1);x=(0,1)+(0,)*(m-2)
 def mul(a,b,h):
  z=[0]*(2*m-1)
  for i in range(m):
   for j in range(m):z[i+j]+=a[i]*b[j]
  for k in range(2*m-2,m-1,-1):
   c=z[k]
   for j in range(m):z[k-m+j]-=c*h[j]
  return tuple(t%4 for t in z[:m])
 def power(a,k,h):
  b=one
  while k:
   if k&1:b=mul(b,a,h)
   a=mul(a,a,h);k//=2
  return b
 lift=None
 for bits in range(1<<m):
  cand=tuple(((base>>i)&1)+2*((bits>>i)&1) for i in range(m))
  if power(x,q-1,cand)==one:lift=cand;break
 if lift is None:raise ValueError('no Teichmuller polynomial lift found')
 T=[zero];v=one
 for _ in range(q-1):T.append(v);v=mul(v,x,lift)
 if v!=one or len(set(T))!=q:raise ValueError('Teichmuller representatives failed')
 def residue(a):return sum((a[i]&1)<<i for i in range(m))
 table={residue(a):a for a in T}
 if len(table)!=q:raise ValueError('Teichmuller residues are not distinct')
 if set(map(tuple,extra))&set(T):raise ValueError('extra evaluations overlap the Teichmuller set')
 def sigma(a):
  u=table[residue(a)];delta=[(a[i]-u[i])%4 for i in range(m)]
  if any(t%2 for t in delta):raise ValueError('invalid Teichmuller decomposition')
  v=table[sum((delta[i]//2)<<i for i in range(m))];u2=mul(u,u,lift);v2=mul(v,v,lift)
  return tuple((u2[i]+2*v2[i])%4 for i in range(m))
 def trace(a):
  out=zero;v=a
  for _ in range(m):out=tuple((out[i]+v[i])%4 for i in range(m));v=sigma(v)
  if v!=a or any(out[1:]):raise ValueError('Galois trace did not land in Z4')
  return out[0]
 basis=[tuple(int(i==j) for i in range(m)) for j in range(m)]
 points=T+list(map(tuple,extra));G=np.array([[trace(mul(a,t,lift)) for t in points] for a in basis]+[[1]*len(points)],np.int64)
 labels=np.arange(count,dtype=np.int64);coeff=np.column_stack([(labels//(4**j))%4 for j in range(m+1)]);words=coeff@G%4;lee=np.array([0,1,2,1],np.int64)[words].sum(1);hist=np.bincount(lee,minlength=2*len(points)+1)
 if hist[0]!=1 or int(hist.sum())!=count:raise ValueError('generator is not injective')
 min_distance=int(np.flatnonzero(hist[1:])[0]+1)
 # The Gray map is an exact Z4 Lee-to-binary Hamming isometry.
 gray=np.array([0,1,3,2],np.int64)
 for a,b in product(range(4),repeat=2):
  if int(gray[a]^gray[b]).bit_count()!=int(np.array([0,1,2,1])[(a-b)%4]):raise ValueError('invalid Gray map')
 data={'m':m,'binary_modulus':base,'lift':list(lift),'teichmuller':[list(t) for t in T],'extra_evaluations':extra,'generator':G.tolist(),'gray_map':gray.tolist(),'lee_histogram':hist.tolist()}
 facts={'quaternary_length':len(points),'binary_length':2*len(points),'cardinality':count,'minimum_distance':min_distance,'ring':'GR(4,%d)'%m,'scope':'exact binary distance via the Z4-linear Lee spectrum'}
 return data,facts,words


def construct(inputs,params,context):
 check_params(params,('m','extra_evaluations'))
 if inputs:raise ValueError('Kerdock construction takes no inputs')
 data,facts,_=_construct(params.get('m'),context,params.get('extra_evaluations'))
 return Artifact('galois_kerdock_code',data,facts,Verification('exact','Teichmuller trace and exhaustive Z4 Lee spectrum',{'facts':facts,'data_sha256':digest(data)}))


def verify_artifact(artifact,context=None):
 context=context or RunContext(Path('.'))
 if artifact.kind!='galois_kerdock_code' or set(artifact.data)!={'m','binary_modulus','lift','teichmuller','extra_evaluations','generator','gray_map','lee_histogram'}:raise ValueError('invalid Kerdock artifact')
 data,facts,_=_construct(artifact.data['m'],context,artifact.data['extra_evaluations'])
 if data!=artifact.data:raise ValueError('Kerdock recipe or spectrum was modified')
 from ..certification import compare_facts
 compare_facts(artifact.facts,facts)
 return Verification('exact','independent replay of Galois-ring trace and exhaustive Lee spectrum',{'facts':facts,'data_sha256':digest(data)})


def spherical(inputs,params,context):
 check_params(params,('coordinates',))
 if len(inputs)!=1 or inputs[0].kind!='galois_kerdock_code':raise ValueError('Kerdock code input required')
 verify_artifact(inputs[0],context);data,facts,words=_construct(inputs[0].data['m'],context,inputs[0].data['extra_evaluations']);length=facts['binary_length'];indices=params.get('coordinates',list(range(length)))
 if not isinstance(indices,list) or not indices or any(type(i) is not int or i<0 or i>=length for i in indices) or len(set(indices))!=len(indices):raise ValueError('coordinates must be distinct binary positions')
 context.limit('max_dimension',len(indices));gray=np.array(data['gray_map'],np.int64)[words];bits=np.empty((len(words),length),np.int8);bits[:,0::2]=(gray>>1)&1;bits[:,1::2]=gray&1;coords=(1-2*bits[:,indices]).astype(int).tolist();out={'coordinates':coords,'norm_squared':len(indices)};report=verify_integer_spherical(out,context)
 return Artifact('spherical_code',out,report.details['facts'],report)


def search_even_extensions(inputs,params,context):
 """Search order-two or arbitrary non-Teichmuller evaluation points."""
 from time import monotonic
 import math
 if len(inputs)!=1 or inputs[0].kind!='galois_kerdock_code':raise ValueError('Kerdock code input required')
 verify_artifact(inputs[0],context);check_params(params,('count','trials','seed','seconds','family'))
 if inputs[0].data['extra_evaluations']:raise ValueError('search requires an unextended source code')
 m=inputs[0].data['m'];q=1<<m;count=params.get('count');trials=params.get('trials',100);seed=params.get('seed',0);seconds=params.get('seconds',5);family=params.get('family','even')
 if family not in ('even','all'):raise ValueError('family must be even or all')
 if type(count) is not int or count<1:raise ValueError('count must be positive')
 if type(trials) is not int or trials<1:raise ValueError('trials must be positive')
 if type(seed) is not int or seed<0:raise ValueError('seed must be nonnegative')
 if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
 context.limit('max_trials',trials);context.limit('solver_seconds',seconds);context.limit('max_dimension',2*(q+count))
 codewords=4**(m+1);base_data,base_facts,base_words=_construct(m,context);lee=np.array([0,1,2,1],np.int16);base=lee[base_words].sum(1).astype(np.int16)
 labels=np.arange(codewords,dtype=np.int64);coeff=np.column_stack([(labels//(4**j))%4 for j in range(m+1)])
 if family=='even':candidates=[[2*((i>>j)&1) for j in range(m)] for i in range(1,q)]
 else:
  teich={tuple(t) for t in base_data['teichmuller']};candidates=[list(v) for v in product(range(4),repeat=m) if v not in teich]
 if count>len(candidates):raise ValueError('count exceeds the evaluation family')
 context.limit('max_rank_work',trials*codewords*2+len(candidates)*codewords*(m+1))
 # Trace linearity computes every extra column from the trace-pairing matrix.
 G=np.asarray(base_data['generator'],np.int64);pairing=G[:m,1:m+1];extra=(pairing@np.asarray(candidates,np.int64).T)%4;extra=np.vstack([extra,np.ones(len(candidates),np.int64)]);C=lee[(coeff@extra)%4].astype(np.int16)
 rng=np.random.default_rng(seed);selected=np.array(rng.choice(len(candidates),count,replace=False),int);weights=base+C[:,selected].sum(1);best=int(weights[1:].min());bestsel=selected.copy();deadline=monotonic()+seconds
 for t in range(trials):
  if monotonic()>=deadline:break
  remove=int(rng.integers(count));available=np.ones(len(candidates),bool);available[selected]=False;options=np.flatnonzero(available)
  if not len(options):break
  add=int(rng.choice(options));nextweights=weights-C[:,selected[remove]]+C[:,add];score=int(nextweights[1:].min())
  if score>=int(weights[1:].min()):selected[remove]=add;weights=nextweights
  if score>best:best=score;bestsel=selected.copy()
 extra=[candidates[i] for i in sorted(bestsel)];data,facts,_=_construct(m,context,extra)
 if facts['minimum_distance']!=best:raise ValueError('search score did not replay')
 return Artifact('galois_kerdock_code',data,facts,Verification('exact','searched evaluation set; exact Galois-ring trace and Lee spectrum',{'facts':facts,'data_sha256':digest(data)}))

def _packed_words(data,coordinates,context):
 rebuilt,facts,words=_construct(data['m'],context,data['extra_evaluations'])
 if rebuilt!=data:raise ValueError('Kerdock source recipe was modified')
 length=facts['binary_length']
 if not isinstance(coordinates,list) or not coordinates or any(type(i) is not int or not 0<=i<length for i in coordinates) or len(set(coordinates))!=len(coordinates):raise ValueError('coordinates must be distinct binary positions')
 d=len(coordinates)
 if d>63:raise ValueError('packed translate backend supports at most 63 selected bits')
 gray=np.array(data['gray_map'],np.uint64)[words];bits=np.empty((len(words),length),np.uint64);bits[:,0::2]=gray&1;bits[:,1::2]=(gray>>1)&1
 masks=np.zeros(len(words),np.uint64)
 for j,i in enumerate(coordinates):masks|=bits[:,i]<<np.uint64(j)
 if len(np.unique(masks))!=len(masks):raise ValueError('puncturing produced duplicate source words')
 return masks,d


def verify_union(artifact,context=None):
 context=context or RunContext(Path('.'));data=artifact.data
 if artifact.kind!='galois_translate_code' or set(data)!={'parent','coordinates','translates','add_d_roots'}:raise ValueError('invalid Galois translate code')
 source=data['parent'];coords=data['coordinates'];trans=data['translates'];roots=data['add_d_roots']
 if not isinstance(source,dict) or type(roots) is not bool:raise ValueError('invalid Galois translate source')
 masks,d=_packed_words(source,coords,context);context.limit('max_dimension',d);N=len(masks);context.limit('max_pairs',N*N)
 if not isinstance(trans,list) or not trans or any(type(t) is not int or t<0 or t>=1<<d for t in trans) or len(set(trans))!=len(trans):raise ValueError('invalid binary translates')
 context.limit('max_codewords',len(trans)*N);D=np.unique((masks[:,None]^masks[None,:]).ravel());nz=D[D!=0];target=(d+3)//4
 if len(nz)==0 or int(np.bitwise_count(nz).min())<target:raise ValueError('punctured source violates kissing distance')
 for i,t in enumerate(trans):
  for u in trans[:i]:
   if int(np.bitwise_count(D^np.uint64(t^u)).min())<target:raise ValueError('translate cross-distance violates kissing condition')
 if roots and d<8:raise ValueError('D-root layer requires at least eight dimensions')
 total=N*len(trans)+(2*d*(d-1) if roots else 0)
 facts={'dimension':d,'ambient_dimension':d,'dimension_scope':'ambient binary embedding; span not asserted','cardinality':total,'binary_codewords':N*len(trans),'translate_layers':len(trans),'minimum_binary_distance_lower_bound':target,'root_points':2*d*(d-1) if roots else 0}
 return Verification('exact','Z4 Lee spectrum, exhaustive binary differences and optional D-root compatibility',{'facts':facts,'data_sha256':digest(data),'difference_count':len(D)})


def find_translates(inputs,params,context):
 from itertools import combinations
 from time import monotonic
 import math
 if len(inputs)!=1 or inputs[0].kind!='galois_kerdock_code':raise ValueError('Kerdock code input required')
 verify_artifact(inputs[0],context);check_params(params,('coordinates','trials','seed','seconds','add_d_roots'))
 data=inputs[0].data;length=2*len(data['generator'][0]);coords=params.get('coordinates',list(range(length)));masks,d=_packed_words(data,coords,context)
 if d>20:raise ValueError('exhaustive translate discovery supports at most 20 bits')
 count=1<<d;context.limit('max_codewords',count);N=len(masks);context.limit('max_pairs',N*N)
 trials=params.get('trials',20);seed=params.get('seed',0);seconds=params.get('seconds',5);roots=params.get('add_d_roots',False)
 if type(trials) is not int or trials<1 or type(seed) is not int or seed<0 or type(roots) is not bool:raise ValueError('invalid translate search parameters')
 if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
 context.limit('max_trials',trials);context.limit('solver_seconds',seconds);target=(d+3)//4;D=np.unique((masks[:,None]^masks[None,:]).ravel());nz=D[D!=0]
 if len(nz)==0 or int(np.bitwise_count(nz).min())<target:raise ValueError('punctured source violates kissing distance')
 errors=[sum(1<<i for i in S) for w in range(target) for S in combinations(range(d),w)];context.limit('max_rank_work',len(errors)*len(D)+trials*count)
 forbidden=np.zeros(count,bool)
 for e in errors:forbidden[(D^np.uint64(e)).astype(np.int64)]=True
 good=~forbidden;allidx=np.arange(count,dtype=np.uint64);rng=np.random.default_rng(seed);best=[0];deadline=monotonic()+seconds;work=len(errors)*len(D)
 for _ in range(trials):
  if monotonic()>=deadline:break
  avail=good.copy();chosen=[0]
  while np.any(avail):
   work+=count;context.limit('max_rank_work',work)
   opts=np.flatnonzero(avail);t=int(rng.choice(opts));chosen.append(t);avail&=good[(allidx^np.uint64(t)).astype(np.int64)]
  if len(chosen)>len(best):best=chosen
 output={'parent':data,'coordinates':coords,'translates':best,'add_d_roots':roots};artifact=Artifact('galois_translate_code',output);report=verify_union(artifact,context);artifact.facts=report.details['facts'];artifact.verification=report;return artifact

def translate_search(inputs,params,context):
 from .binary_translate_search import binary_translate_search
 if len(inputs)!=1 or inputs[0].kind!='galois_kerdock_code':raise ValueError('Kerdock code input required')
 verify_artifact(inputs[0],context)
 allowed=('alphabet','references','reference_offsets','candidate_offset','minimum_distance','rounds','oracle_steps','exact_pair_limit','seed','seconds')
 check_params(params,allowed);p=dict(params);p.update(generator=inputs[0].data['generator'],modulus=4,block_map=inputs[0].data['gray_map'],bit_width=2)
 return binary_translate_search([],p,context)


OPERATORS={'galois.kerdock':construct,'galois.spherical':spherical,'galois.translate_search':translate_search,'galois.evaluation_search':search_even_extensions,'galois.find_translates':find_translates}
INPUT_KINDS={'galois.kerdock':(),'galois.spherical':('galois_kerdock_code',),'galois.translate_search':('galois_kerdock_code',),'galois.evaluation_search':('galois_kerdock_code',),'galois.find_translates':('galois_kerdock_code',)}
OUTPUT_KINDS={'galois.kerdock':('galois_kerdock_code',),'galois.spherical':('spherical_code',),'galois.translate_search':('search_result',),'galois.evaluation_search':('galois_kerdock_code',),'galois.find_translates':('galois_translate_code',)}
