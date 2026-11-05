"""Implicit binary codes from affine subspace cosets, exactly checked by span enumeration."""
from ..core import Artifact,RunContext,Verification,check_params,digest,integer


def verify_artifact(a,context):
 if a.kind!='binary_code' or set(a.data)!={'backend','length','generators','representatives'} or a.data['backend']!='binary_span_cosets':raise ValueError('invalid span-coset code')
 p=a.data;n=integer(p['length'],'length');context.limit('max_dimension',n);G=p['generators'];R=p['representatives']
 if not isinstance(G,list) or not isinstance(R,list) or not R or any(type(v) is not int or v<0 or v.bit_length()>n for v in G+R):raise ValueError('invalid binary vectors')
 k=len(G)
 if k>30:raise ValueError('span rank exceeds hard bound')
 count=1<<k;context.limit('max_codewords',count);context.limit('max_pairs',count*(1+len(R)*(len(R)-1)//2));context.limit('max_rank_work',count*max(1,len(R)*(len(R)-1)//2))
 H=[0]
 for g in G:H += [v^g for v in H]
 if len(set(H))!=count:raise ValueError('dependent span generators')
 kernel=min((v.bit_count() for v in H if v),default=n+1);distance=kernel
 for i,a0 in enumerate(R):
  for b0 in R[:i]:
   delta=a0^b0;d=min((v^delta).bit_count() for v in H)
   if d==0:raise ValueError('duplicate affine cosets')
   distance=min(distance,d)
 facts={'length':n,'cardinality':count*len(R),'min_distance':distance if count*len(R)>1 else None}
 from ..certification import compare_facts
 compare_facts(a.facts,facts)
 return Verification('exact','all span words and inter-coset differences enumerated exactly',{'facts':facts,'span_size':count,'coset_count':len(R),'data_sha256':digest(p)})


def build(inputs,params,context):
 check_params(params,('length','generators','representatives'))
 if inputs:raise ValueError('binary span cosets takes no inputs')
 data={'backend':'binary_span_cosets','length':params.get('length'),'generators':params.get('generators'),'representatives':params.get('representatives')};a=Artifact('binary_code',data);a.verification=verify_artifact(a,context);a.facts=a.verification.details['facts'];return a

OPERATORS={'codes.binary_span_cosets':build}
INPUT_KINDS={'codes.binary_span_cosets':()}
OUTPUT_KINDS={'codes.binary_span_cosets':('binary_code',)}
