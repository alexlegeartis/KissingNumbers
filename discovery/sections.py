"""Exact sections inheriting separation from any verified coordinate backend."""
from fractions import Fraction
from math import lcm

import numpy as np

from .core import Artifact,Verification,integer,digest
from .exact_gram import bounded_fraction,check_fraction_size,exact_integer_rank,rational_json


def _depth(source,context):
    count=0;current=source
    while isinstance(current,dict) and isinstance(current.get('data'),dict):
        count+=1;context.limit('max_nodes',count)
        data=current['data']
        if data.get('backend')!='certified_coordinate_section':break
        current=data.get('source')


def _select(rows,weights,normals,context):
    n,d=len(rows),len(rows[0])
    if not isinstance(normals,list) or not normals:
        raise ValueError('Nonempty list of section normals required')
    context.limit('max_nodes',len(normals))
    context.limit('max_gram_entries',len(normals)*d)
    context.limit('max_rank_work',3*n*d*len(normals))
    integral=[];canonical=[]
    for raw in normals:
        if not isinstance(raw,list) or len(raw)!=d:raise ValueError('Normal width mismatch')
        values=[bounded_fraction(x,context,'section normal') for x in raw]
        if not any(values):raise ValueError('Section normal must be nonzero')
        denominator=1
        for value in values:
            denominator=lcm(denominator,value.denominator)
            context.limit('max_rational_bits',denominator.bit_length())
        coefficients=[int(value*denominator)*weight for value,weight in zip(values,weights)]
        for value in coefficients:context.limit('max_rational_bits',abs(value).bit_length())
        integral.append(coefficients);canonical.append([rational_json(x) for x in values])
    largest=max(abs(x) for row in rows for x in row)
    bound=largest*max(sum(abs(x) for x in normal) for normal in integral)
    dtype=np.int64 if bound<2**62 else object
    matrix=np.asarray(integral,dtype=dtype)
    selected=[]
    for at in range(0,n,128):
        dots=np.asarray(rows[at:at+128],dtype=dtype)@matrix.T
        selected.extend((np.flatnonzero(np.all(dots==0,axis=1))+at).tolist())
    if not selected:raise ValueError('Section contains no points')
    normal_rank,_=exact_integer_rank(integral,context)
    return selected,canonical,normal_rank


def verify_certified_section(artifact,context):
    from .certification import verify_artifact,compare_facts
    data=artifact.data
    if artifact.kind!='spherical_code' or not isinstance(data,dict) or set(data)!={'backend','coordinates','weights','norm_squared','source','normals'} or data['backend']!='certified_coordinate_section':
        raise ValueError('Unknown certified section descriptor')
    _depth(data['source'],context)
    source=Artifact.from_dict(data['source'])
    rows=source.data.get('coordinates')
    if source.kind!='spherical_code' or not isinstance(rows,list) or not rows or not isinstance(rows[0],list) or not rows[0]:
        raise ValueError('Section source needs verified explicit coordinates')
    n,d=len(rows),len(rows[0]);context.limit('max_materialized_points',n);context.limit('max_dimension',d)
    if not isinstance(data['coordinates'],list):raise ValueError('Section output coordinates must be a list')
    context.limit('max_materialized_points',len(data['coordinates']))
    for row in data['coordinates']:
        if not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row):
            raise ValueError('Section output coordinates must be exact integer rows')
    if not isinstance(data['weights'],list) or len(data['weights'])!=d:
        raise ValueError('Section output metric width')
    for weight in data['weights']:integer(weight,'positive output metric weight')
    integer(data['norm_squared'],'positive output norm')
    weights=source.data.get('weights',[1]*d)
    if not isinstance(weights,list) or len(weights)!=d:raise ValueError('Section source metric width')
    for weight in weights:integer(weight,'positive section metric weight')
    # Shape/normal guards before recursive source proof work.
    normals=data['normals']
    if not isinstance(normals,list) or not normals or any(not isinstance(row,list) or len(row)!=d for row in normals):
        raise ValueError('Nonempty correctly shaped normals required')
    context.limit('max_nodes',len(normals))
    for row in normals:
        values=[bounded_fraction(x,context,'section normal') for x in row]
        if not any(values):raise ValueError('Section normal must be nonzero')
    proof=verify_artifact(source,context)
    if proof.status!='exact':raise ValueError('Section source must independently verify exactly')
    source_facts=proof.details.get('facts',proof.details)
    selected,canonical,normal_rank=_select(rows,weights,normals,context)
    expected=[rows[i] for i in selected]
    if data['coordinates']!=expected or data['weights']!=weights or data['norm_squared']!=source.data['norm_squared']:
        raise ValueError('Section coordinates, norm or metric differ from the complete source subset')
    source_rank=integer(source_facts.get('rank_exact'),'exact source rank')
    upper=min(source_rank,d-normal_rank)
    lower,work=exact_integer_rank(expected[:min(len(expected),max(64,4*d))],context)
    if lower<upper:lower,work=exact_integer_rank(expected,context)
    if lower>upper:raise ValueError('Section rank exceeds exact normal/source bound')
    facts={'dimension':lower,'rank_exact':lower,'ambient_dimension':d,'cardinality':len(expected),
           'norm_squared':data['norm_squared']}
    compare_facts(artifact.facts,facts)
    return Verification('exact','complete exact central section membership; separation inherited from independently verified source',
        {'facts':facts,'source_cardinality':n,'source_data_sha256':digest(source.data),
         'source_verification_method':proof.method,'normal_rank':normal_rank,
         'rank_scalar_work':work,'source_separation_inherited':True,
         'unordered_output_pairs_covered':len(expected)*(len(expected)-1)//2,
         'normals':canonical,'selected_indices':selected,'data_sha256':digest(data)})


def certified_section(source,params,context):
    if ('normal' in params)==('normals' in params):raise ValueError('Provide exactly one of normal or normals')
    rows=source.data.get('coordinates')
    if not isinstance(rows,list) or not rows or not isinstance(rows[0],list) or not rows[0]:raise ValueError('Explicit section source required')
    context.limit('max_materialized_points',len(rows));context.limit('max_dimension',len(rows[0]))
    _depth({'data':source.data},context)
    # Geometry validation belongs to the independent source checker. Selection
    # is repeated after that replay by the terminal checker below.
    normals=params.get('normals',[params.get('normal')])
    if not isinstance(normals,list) or not normals or any(not isinstance(row,list) or len(row)!=len(rows[0]) for row in normals):
        raise ValueError('Nonempty correctly shaped normals required')
    context.limit('max_nodes',len(normals))
    for row in normals:
        values=[bounded_fraction(x,context,'section normal') for x in row]
        if not any(values):raise ValueError('Section normal must be nonzero')
    from .certification import verify_artifact
    if verify_artifact(source,context).status!='exact':raise ValueError('Exactly verified section source required')
    weights=source.data.get('weights',[1]*len(rows[0]))
    indices,canonical,_=_select(rows,weights,normals,context)
    data={'backend':'certified_coordinate_section','coordinates':[rows[i][:] for i in indices],
          'weights':weights[:],'norm_squared':source.data['norm_squared'],
          'source':source.to_dict(),'normals':canonical}
    result=Artifact('spherical_code',data)
    result.verification=verify_certified_section(result,context);result.facts=result.verification.details['facts']
    return result
