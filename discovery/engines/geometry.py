"""Exact transfer of verified binary distance to implicit spherical geometry."""
from ..core import Artifact, Verification, check_params, integer


def verify_sign_embedding(artifact,context):
    from ..certification import verify_artifact,compare_facts
    data=artifact.data
    if set(data)!= {'embedding','binary_code'} or data['embedding']!='binary_sign':
        raise ValueError('Unknown implicit embedding recipe')
    source=Artifact.from_dict(data['binary_code'])
    if source.kind!='binary_code': raise ValueError('Sign embedding requires a binary code')
    report=verify_artifact(source,context)
    if report.status!='exact': raise ValueError('Binary distance is not exactly certified')
    facts=report.details.get('facts',report.details)
    N=integer(facts.get('length'),'binary length');M=integer(facts.get('cardinality'),'code cardinality')
    context.limit('max_dimension',N)
    distance=integer(facts.get('min_distance'),'binary distance',minimum=0)
    if M>1 and distance==0: raise ValueError('Distinctness requires positive minimum distance')
    if 4*distance<N: raise ValueError('Binary distance does not imply the kissing inequality')
    # N is an ambient dimension, not a computed exact rank. An upper bound on
    # dimension suffices for a construction in R^N; do not claim full span.
    derived={'dimension':N,'ambient_dimension':N,'cardinality':M,'norm_squared':N,
             'max_inner_product_upper_bound':N-2*distance,
             'dimension_scope':'ambient; exact span rank not computed'}
    compare_facts(artifact.facts,derived)
    return Verification('exact','binary sign embedding from independently reverified distance',
                        {'facts':derived,'binary_verification':report.to_dict(),'materialized':False})


def sign_embed(inputs,params,context):
    check_params(params,())
    artifact=Artifact('implicit_spherical_code',{'embedding':'binary_sign','binary_code':inputs[0].to_dict()})
    artifact.verification=verify_sign_embedding(artifact,context)
    artifact.facts=artifact.verification.details['facts']
    return artifact


def integer_code(inputs,params,context):
    from ..certification import verify_integer_spherical
    check_params(params,('coordinates','weights','norm_squared'))
    data=dict(params)
    artifact=Artifact('spherical_code',data)
    artifact.verification=verify_integer_spherical(data,context)
    artifact.facts=artifact.verification.details['facts']
    return artifact

def orthogonal_union(inputs,params,context):
    from ..certification import verify_artifact,verify_integer_spherical
    check_params(params,())
    if any(a.kind!='spherical_code' or 'coordinates' not in a.data for a in inputs):
        raise ValueError('orthogonal_union currently requires two explicit integer spherical codes')
    for source in inputs:
        if verify_artifact(source,context).status!='exact':raise ValueError('Exact inputs required')
    a,b=(x.data for x in inputs)
    da,db=len(a['coordinates'][0]),len(b['coordinates'][0])
    context.limit('max_materialized_points',len(a['coordinates'])+len(b['coordinates']))
    context.limit('max_dimension',da+db)
    n1,n2=a['norm_squared'],b['norm_squared']
    coords=[row+[0]*db for row in a['coordinates']]+[[0]*da+row for row in b['coordinates']]
    weights=[w*n2 for w in a.get('weights',[1]*da)]+[w*n1 for w in b.get('weights',[1]*db)]
    result=Artifact('spherical_code',{'coordinates':coords,'weights':weights,'norm_squared':n1*n2})
    result.verification=verify_integer_spherical(result.data,context)
    result.facts=result.verification.details['facts']
    return result

OPERATORS={'geometry.sign_embed':sign_embed,'geometry.integer_code':integer_code,
           'geometry.orthogonal_union':orthogonal_union}
INPUT_KINDS={'geometry.sign_embed':('binary_code',),'geometry.integer_code':(),
             'geometry.orthogonal_union':('spherical_code','spherical_code')}
