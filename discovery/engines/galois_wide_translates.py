"""Random binary translates of nonlinear Galois trace codes, with exact replay.

The Gray image is nonlinear: cross-distances are checked against the complete
binary XOR difference set, never against quaternary message differences alone.
"""
from pathlib import Path
from time import monotonic
import math,numpy as np
from ..core import Artifact,RunContext,Verification,check_params,digest
from .galois_kerdock import _construct,verify_artifact as verify_parent


def _difference_set(data,context):
 parent=Artifact('galois_kerdock_code',data);verify_parent(parent,context);rebuilt,facts,words=_construct(data['m'],context,data['extra_evaluations']);d=facts['binary_length'];N=len(words)
 if d>128:raise ValueError('wide translate backend supports at most 128 bits')
 context.limit('max_pairs',N*N);context.limit('max_gram_entries',2*N*N);context.limit('max_rank_work',2*N*N)
 gray=np.array(data['gray_map'],np.uint64)[words];lo=np.zeros(N,np.uint64);hi=np.zeros(N,np.uint64)
 for j in range(gray.shape[1]):
  if j<32:lo|=gray[:,j]<<np.uint64(2*j)
  else:hi|=gray[:,j]<<np.uint64(2*(j-32))
 pairs=np.empty(N*N,dtype=[('lo','u8'),('hi','u8')]);pairs['lo']=(lo[:,None]^lo[None,:]).ravel();pairs['hi']=(hi[:,None]^hi[None,:]).ravel();D=np.unique(pairs);context.limit('max_codewords',len(D))
 return D,d,facts


def _valid_shift(t,d):
 if not isinstance(t,list) or len(t)!=2 or any(type(x) is not int or x<0 for x in t):return False
 if t[0]>=1<<min(d,64):return False
 return t[1]==0 if d<=64 else t[1]<(1<<(d-64))


def verify_artifact(artifact,context=None):
 context=context or RunContext(Path('.'));data=artifact.data
 if artifact.kind!='galois_wide_translate_code' or set(data)!={'parent','shifts','add_d_roots'}:raise ValueError('invalid wide translate artifact')
 if not isinstance(data['parent'],dict) or type(data['add_d_roots']) is not bool:raise ValueError('invalid wide source')
 D,d,parentfacts=_difference_set(data['parent'],context);shifts=data['shifts'];target=(d+3)//4
 if parentfacts['minimum_distance']<target:raise ValueError('parent code violates kissing distance')
 if not isinstance(shifts,list) or not shifts or any(not _valid_shift(t,d) for t in shifts) or len(set(map(tuple,shifts)))!=len(shifts):raise ValueError('invalid binary shifts')
 L=len(shifts);context.limit('max_codewords',L*parentfacts['cardinality']);context.limit('max_rank_work',max(2*parentfacts['cardinality']**2,L*(L-1)//2*len(D)))
 lo,hi=D['lo'],D['hi']
 for i,(a,b) in enumerate(shifts):
  for c,e in shifts[:i]:
   if int((np.bitwise_count(lo^np.uint64(a^c))+np.bitwise_count(hi^np.uint64(b^e))).min())<target:raise ValueError('wide translate cross-distance violates kissing condition')
 if data['add_d_roots'] and d<8:raise ValueError('D-root layer requires at least eight dimensions')
 roots=2*d*(d-1) if data['add_d_roots'] else 0;facts={'dimension':d,'ambient_dimension':d,'dimension_scope':'ambient binary embedding; span not asserted','cardinality':L*parentfacts['cardinality']+roots,'binary_codewords':L*parentfacts['cardinality'],'translate_layers':L,'minimum_binary_distance_lower_bound':target,'root_points':roots}
 return Verification('exact','exhaustive nonlinear binary difference set and exact translate cross-distances',{'facts':facts,'data_sha256':digest(data),'difference_count':len(D)})


def search(inputs,params,context):
 if len(inputs)!=1 or inputs[0].kind!='galois_kerdock_code':raise ValueError('Kerdock code input required')
 check_params(params,('trials','seed','seconds','max_layers','add_d_roots'))
 trials=params.get('trials',200);seed=params.get('seed',0);seconds=params.get('seconds',10);max_layers=params.get('max_layers',100);roots=params.get('add_d_roots',False)
 if type(trials) is not int or trials<1 or type(seed) is not int or seed<0 or type(max_layers) is not int or max_layers<1 or type(roots) is not bool:raise ValueError('invalid wide translate search parameters')
 if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:raise ValueError('seconds must be positive')
 context.limit('max_trials',trials);context.limit('solver_seconds',seconds);D,d,facts=_difference_set(inputs[0].data,context);target=(d+3)//4
 if facts['minimum_distance']<target:raise ValueError('parent code violates kissing distance')
 if roots and d<8:raise ValueError('D-root layer requires at least eight dimensions')
 rng=np.random.default_rng(seed);shifts=[[0,0]];deadline=monotonic()+seconds;work=2*facts['cardinality']**2;lo,hi=D['lo'],D['hi']
 exhausted=False
 for _ in range(trials):
  if monotonic()>=deadline or len(shifts)>=max_layers or (len(shifts)+1)*facts['cardinality']>context.budget['max_codewords']:break
  a=int(rng.integers(0,1<<min(d,64),dtype=np.uint64));b=0 if d<=64 else int(rng.integers(0,1<<(d-64),dtype=np.uint64));ok=True
  for c,e in shifts:
   if work+len(D)>context.budget['max_rank_work']:exhausted=True;ok=False;break
   work+=len(D)
   if int((np.bitwise_count(lo^np.uint64(a^c))+np.bitwise_count(hi^np.uint64(b^e))).min())<target:ok=False;break
  if exhausted:break
  if ok:shifts.append([a,b])
 output={'parent':inputs[0].data,'shifts':shifts,'add_d_roots':roots};artifact=Artifact('galois_wide_translate_code',output);report=verify_artifact(artifact,context);artifact.facts=report.details['facts'];artifact.verification=report;return artifact


OPERATORS={'galois.random_translates':search}
INPUT_KINDS={'galois.random_translates':('galois_kerdock_code',)}
OUTPUT_KINDS={'galois.random_translates':('galois_wide_translate_code',)}
