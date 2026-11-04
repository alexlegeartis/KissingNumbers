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


def gram_code(inputs,params,context):
    """Exact symbolic real spherical code from an arbitrary rational Gram."""
    from ..exact_gram import verify_rational_gram
    check_params(params,('gram','norm_squared'))
    if inputs: raise ValueError('geometry.gram_code takes no inputs')
    artifact=Artifact('gram_spherical_code',dict(params))
    artifact.verification=verify_rational_gram(artifact.data,context)
    artifact.facts=artifact.verification.details['facts']
    return artifact


def _section_rows(source,normal,context):
    from ..exact_gram import bounded_fraction,check_fraction_size,rational_json
    rows=source['coordinates'];d=len(rows[0])
    if not isinstance(normal,list) or len(normal)!=d:
        raise ValueError('Section normal must match the ambient coordinate width')
    values=[bounded_fraction(x,context,'section normal') for x in normal]
    if not any(values):raise ValueError('Section normal must be nonzero')
    weights=source.get('weights',[1]*d)
    coefficients=[check_fraction_size(value*weight,context) for value,weight in zip(values,weights)]
    selected=[]
    for i,row in enumerate(rows):
        dot=0
        for coefficient,value in zip(coefficients,row):
            dot=check_fraction_size(dot+coefficient*value,context)
        if dot==0:selected.append(i)
    if not selected:raise ValueError('Hyperplane section contains no points')
    return selected,[rational_json(x) for x in values]


def verify_section(artifact,context):
    """Replay exact source geometry and the full weighted hyperplane predicate."""
    from ..certification import verify_integer_spherical,compare_facts
    from ..core import digest
    data=artifact.data
    if set(data)!={'coordinates','weights','norm_squared','section'}:
        raise ValueError('Unknown explicit section descriptor fields')
    descriptor=data['section']
    if not isinstance(descriptor,dict) or set(descriptor)!={'source','normal'}:
        raise ValueError('Section needs exact source and normal')
    source=descriptor['source']
    if not isinstance(source,dict) or set(source)!={'coordinates','weights','norm_squared'}:
        raise ValueError('Section source must be a plain exact weighted coordinate code')
    # Reject the normal shape before any source pair/rank work.
    rows=source.get('coordinates')
    if not isinstance(rows,list) or not rows or not isinstance(rows[0],list):
        raise ValueError('Section source needs nonempty explicit coordinates')
    context.limit('max_materialized_points',len(rows))
    context.limit('max_dimension',len(rows[0]))
    context.limit('max_pairs',len(rows)*(len(rows)-1)//2)
    normal=descriptor['normal']
    if not isinstance(normal,list) or len(normal)!=len(rows[0]):
        raise ValueError('Section normal must match the ambient coordinate width')
    from ..exact_gram import bounded_fraction
    values=[bounded_fraction(x,context,'section normal') for x in normal]
    if not any(values):raise ValueError('Section normal must be nonzero')
    source_report=verify_integer_spherical(source,context)
    indices,canonical_normal=_section_rows(source,normal,context)
    expected=[source['coordinates'][i] for i in indices]
    if data['coordinates']!=expected or data['weights']!=source['weights'] or data['norm_squared']!=source['norm_squared']:
        raise ValueError('Stored section differs from exact selected source rows or metric')
    report=verify_integer_spherical(data,context)
    compare_facts(artifact.facts,report.details['facts'])
    report.method='exact rational weighted central hyperplane selection and independent integer geometry'
    report.details['section']={'normal':canonical_normal,'normal_convention':'sum(weights[j]*normal[j]*coordinate[j]) == 0',
                               'selected_indices':indices,'source_data_sha256':digest(source),
                               'source_rank_exact':source_report.details['facts']['rank_exact'],
                               'source_cardinality':source_report.details['facts']['cardinality'],
                               'ambient_metric_preserved':True}
    return report


def section(inputs,params,context):
    from ..certification import verify_artifact
    check_params(params,('normal',))
    if len(inputs)!=1 or inputs[0].kind!='spherical_code' or 'coordinates' not in inputs[0].data:
        raise ValueError('geometry.section requires one explicit exact integer coordinate code')
    source=inputs[0]
    normal=params.get('normal');rows=source.data.get('coordinates')
    if not isinstance(rows,list) or not rows or not isinstance(rows[0],list):
        raise ValueError('Section source needs explicit coordinate rows')
    context.limit('max_materialized_points',len(rows))
    context.limit('max_dimension',len(rows[0]))
    context.limit('max_pairs',len(rows)*(len(rows)-1)//2)
    if not isinstance(normal,list) or len(normal)!=len(rows[0]):
        raise ValueError('Section normal must match the ambient coordinate width')
    from ..exact_gram import bounded_fraction
    checked_normal=[bounded_fraction(x,context,'section normal') for x in normal]
    if not any(checked_normal):raise ValueError('Section normal must be nonzero')
    if verify_artifact(source,context).status!='exact':raise ValueError('Section needs an exactly verified source')
    plain={'coordinates':[row[:] for row in rows],
           'weights':source.data.get('weights',[1]*len(rows[0]))[:],
           'norm_squared':source.data['norm_squared']}
    indices,canonical_normal=_section_rows(plain,normal,context)
    result=Artifact('spherical_code',{'coordinates':[plain['coordinates'][i][:] for i in indices],
                                    'weights':plain['weights'][:],'norm_squared':plain['norm_squared'],
                                    'section':{'source':plain,'normal':canonical_normal}})
    result.verification=verify_section(result,context)
    result.facts=result.verification.details['facts']
    return result

OPERATORS={'geometry.sign_embed':sign_embed,'geometry.integer_code':integer_code,
           'geometry.orthogonal_union':orthogonal_union,'geometry.gram_code':gram_code,
           'geometry.section':section}
INPUT_KINDS={'geometry.sign_embed':('binary_code',),'geometry.integer_code':(),
             'geometry.orthogonal_union':('spherical_code','spherical_code'),'geometry.gram_code':(),
             'geometry.section':('spherical_code',)}
