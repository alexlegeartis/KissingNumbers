"""Exact rational orthogonal orbits with a compressed separation certificate.

Closure, transitivity and one seed against every orbit point cover every pair.
No irreducibility or trivial-stabilizer assumption is made.
"""
from math import isqrt
from pathlib import Path

import numpy as np

from ..core import Artifact,Verification,RunContext,check_params,digest,integer
from ..exact_gram import exact_integer_rank


def _generators(raw,d,context):
    if not isinstance(raw,list) or not raw:
        raise ValueError('Exact orbit generators must be nonempty')
    context.limit('max_nodes',len(raw))
    context.limit('max_gram_entries',len(raw)*d*d)
    context.limit('max_rank_work',3*len(raw)*d*d*d)
    result=[]
    for item in raw:
        if not isinstance(item,dict) or set(item)!={'matrix','denominator'}:
            raise ValueError('Generator needs integer matrix and positive denominator')
        matrix=item['matrix'];den=integer(item['denominator'],'generator denominator')
        context.limit('max_rational_bits',den.bit_length())
        if not isinstance(matrix,list) or len(matrix)!=d or any(not isinstance(row,list) or len(row)!=d for row in matrix):
            raise ValueError('Exact generator dimension mismatch')
        for row in matrix:
            if any(type(x) is not int for x in row):raise ValueError('Generator coefficients must be exact integers')
            for x in row:context.limit('max_rational_bits',abs(x).bit_length())
        for i in range(d):
            for j in range(i+1):
                value=sum(matrix[k][i]*matrix[k][j] for k in range(d))
                context.limit('max_rational_bits',abs(value).bit_length())
                if value!=(den*den if i==j else 0):raise ValueError('Generator is not exactly orthogonal')
        result.append((matrix,den))
    return result


def _image(row,matrix,den,context,sparse=None):
    values=([sum(coefficient*row[j] for j,coefficient in line) for line in sparse]
            if sparse is not None else [sum(x*y for x,y in zip(line,row)) for line in matrix])
    for x in values:
        context.limit('max_rational_bits',abs(x).bit_length())
        if x%den:raise ValueError('Image leaves the supplied integral coordinate scale; not an orbit impossibility proof')
    return tuple(x//den for x in values)


def integer_orbit(inputs,params,context):
    check_params(params,('seed','generators'))
    if inputs:raise ValueError('groups.integer_orbit takes exact generator and seed parameters')
    seed=params.get('seed')
    if not isinstance(seed,list) or not seed or any(type(x) is not int for x in seed):
        raise ValueError('Exact nonzero integer seed required')
    d=len(seed);context.limit('max_dimension',d)
    for x in seed:context.limit('max_rational_bits',abs(x).bit_length())
    norm=sum(x*x for x in seed)
    if not norm:raise ValueError('Exact nonzero integer seed required')
    context.limit('max_rational_bits',norm.bit_length())
    generators=_generators(params.get('generators'),d,context)
    sparse=[[[ (j,x) for j,x in enumerate(line) if x] for line in matrix] for matrix,_ in generators]
    rows=[tuple(seed)];seen={tuple(seed)}
    for row in rows:
        context.limit('max_rank_work',len(rows)*len(generators)*d*d)
        for (matrix,den),nonzero in zip(generators,sparse):
            image=_image(row,matrix,den,context,nonzero)
            if image not in seen:
                # One bad distinct image excludes this full generated orbit;
                # the exception is not a rejection of partial orbit selections.
                if 2*sum(x*y for x,y in zip(seed,image))>norm:
                    raise ValueError('Full generated orbit has an exact distinct-image kissing violation; partial orbits remain open')
                context.limit('max_materialized_points',len(rows)+1)
                context.limit('max_orbit_points',len(rows)+1)
                context.limit('max_pairs',len(rows))
                seen.add(image);rows.append(image)
    artifact=Artifact('spherical_code',{'backend':'exact_integer_orbit','coordinates':[list(row) for row in rows],
        'norm_squared':norm,'orbit_generators':params['generators'],'seed_index':0})
    artifact.verification=verify_artifact(artifact,context)
    artifact.facts=artifact.verification.details['facts']
    return artifact


def verify_artifact(artifact,context=None):
    """Independent block action replay, graph connectivity and seed separation."""
    context=context or RunContext(Path('.'));data=artifact.data
    if artifact.kind!='spherical_code' or not isinstance(data,dict) or set(data)!={'backend','coordinates','norm_squared','orbit_generators','seed_index'} or data['backend']!='exact_integer_orbit':
        raise ValueError('Unknown exact integer orbit descriptor')
    rows=data['coordinates'];norm=integer(data['norm_squared'],'orbit norm')
    if not isinstance(rows,list) or not rows or not isinstance(rows[0],list) or not rows[0]:
        raise ValueError('Nonempty exact orbit rows required')
    n,d=len(rows),len(rows[0]);context.limit('max_dimension',d)
    context.limit('max_materialized_points',n);context.limit('max_orbit_points',n)
    context.limit('max_pairs',n-1)
    seed_index=integer(data['seed_index'],'seed index',minimum=0)
    if seed_index>=n:raise ValueError('Seed index out of range')
    context.limit('max_rational_bits',norm.bit_length())
    for row in rows:
        if not isinstance(row,list) or len(row)!=d or any(type(x) is not int for x in row):raise ValueError('Exact rectangular orbit rows required')
        for x in row:context.limit('max_rational_bits',abs(x).bit_length())
        if sum(x*x for x in row)!=norm:raise ValueError('Orbit norm mismatch')
    lookup={tuple(row):i for i,row in enumerate(rows)}
    if len(lookup)!=n:raise ValueError('Duplicate orbit points')
    generators=_generators(data['orbit_generators'],d,context)
    context.limit('max_rank_work',n*len(generators)*d*d)
    context.limit('max_gram_entries',n*len(generators))
    context.limit('max_gram_entries',min(128,n)*d)
    edges=np.empty((n,len(generators)),dtype=np.int64)
    for column,(matrix,den) in enumerate(generators):
        bound=isqrt(norm)*max(sum(abs(x) for x in line) for line in matrix)
        dtype=np.int64 if bound<2**62 and den<2**62 else object
        value=np.asarray(matrix,dtype=dtype)
        for at in range(0,n,128):
            numerator=np.asarray(rows[at:at+128],dtype=dtype)@value.T
            for local,row in enumerate(numerator):
                if any(int(x)%den for x in row):raise ValueError('Nonintegral image in stored orbit')
                key=tuple(int(x)//den for x in row)
                if key not in lookup:raise ValueError('Stored orbit is not closed')
                edges[at+local,column]=lookup[key]
    visited={seed_index};queue=[seed_index]
    for i in queue:
        for value in edges[i]:
            j=int(value)
            if j not in visited:visited.add(j);queue.append(j)
    if len(visited)!=n:raise ValueError('Stored point set is not one transitive orbit')
    seed=rows[seed_index];maximum=None
    for i,row in enumerate(rows):
        if i==seed_index:continue
        dot=sum(x*y for x,y in zip(seed,row))
        if 2*dot>norm:raise ValueError('Seed has an exact distinct orbit conflict')
        maximum=dot if maximum is None else max(maximum,dot)
    rank,work=exact_integer_rank(rows[:min(n,max(64,4*d))],context)
    if rank<d:
        rank,work=exact_integer_rank(rows,context)
    facts={'dimension':rank,'rank_exact':rank,'ambient_dimension':d,'cardinality':n,
           'norm_squared':norm,'max_inner_product':maximum,'finite_group_proved':rank==d}
    from ..certification import compare_facts
    compare_facts(artifact.facts,facts)
    return Verification('exact','exact orthogonal action closure, transitivity, seed separation and rational rank',
        {'facts':facts,'seed_pairs_checked':n-1,'unordered_pairs_covered':n*(n-1)//2,
         'generator_images_checked':n*len(generators),'rank_scalar_work':work,'data_sha256':digest(data),
         'proof':'Any pair can be moved to a seed pair by an inverse word in the exactly orthogonal generators. Closure and reachability are checked explicitly.',
         'finite_group_scope':'The generated matrix group is finite if the orbit spans the whole ambient space; no abstract group order or irreducibility is claimed.'})


OPERATORS={'groups.integer_orbit':integer_orbit}
INPUT_KINDS={'groups.integer_orbit':()}
OUTPUT_KINDS={'groups.integer_orbit':('spherical_code',)}
