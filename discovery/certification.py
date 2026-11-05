"""Exact coordinate verification and proof dispatch. No trust in saved status."""
from pathlib import Path
import numpy as np
from .core import Artifact, Verification, RunContext, integer, digest, canonical_json


def verify_integer_spherical(data,context=None):
    context=context or RunContext(Path('.'))
    if not isinstance(data,dict): raise ValueError('Coordinates require an object')
    rows=data.get('coordinates'); weights=data.get('weights'); norm=data.get('norm_squared')
    integer(norm,'norm_squared')
    if not isinstance(rows,list) or not rows or not isinstance(rows[0],list) or not rows[0]:
        raise ValueError('Nonempty integer coordinates required')
    n,d=len(rows),len(rows[0]); context.limit('max_materialized_points',n);context.limit('max_dimension',d)
    context.limit('max_pairs',n*(n-1)//2)
    if weights is None: weights=[1]*d
    if not isinstance(weights,list) or len(weights)!=d:
        raise ValueError('Wrong diagonal quadratic form')
    for w in weights: integer(w,'quadratic-form weight')
    for row in rows:
        if not isinstance(row,list) or len(row)!=d:
            raise ValueError('Ragged coordinates')
        if any(isinstance(x,bool) or not isinstance(x,int) for x in row):
            raise ValueError('Exact coordinates must be integers, not rounded floats')
    if len(set(map(tuple,rows)))!=n: raise ValueError('Duplicate points')
    for value in [norm]+weights:
        context.limit('max_rational_bits',value.bit_length())
    for row in rows:
        for value in row:context.limit('max_rational_bits',value.bit_length())
    maxabs=max(abs(x) for row in rows for x in row)
    safe=maxabs*maxabs*sum(weights)<2**62
    dtype=np.int64 if safe else object
    C=np.array(rows,dtype=dtype);W=np.array(weights,dtype=dtype)
    norms=np.sum(C*C*W,axis=1)
    if any(int(x)!=norm for x in norms): raise ValueError('Unequal or incorrect norms')
    context.limit('max_gram_entries',min(128,n)*n)
    # Every dot-product operand, product and partial sum is an integer whose
    # absolute value is bounded by maxabs**2*sum(weights). Below 2**24 these
    # sums are represented exactly in binary32, in any accumulation order.
    binary32=maxabs*maxabs*sum(weights)<2**24
    products=C.astype(np.float32) if binary32 else C
    product_weights=W.astype(np.float32) if binary32 else W
    maximum=None
    for start in range(0,n,128):
        G=(products[start:start+128]*product_weights)@products.T
        for k,row in enumerate(G):
            off=row[start+k+1:]
            if len(off):
                value=int(np.max(off));maximum=value if maximum is None else max(maximum,value)
                if 2*value>norm: raise ValueError('Pairwise kissing inequality violated')
    from .exact_gram import exact_integer_rank
    rank,rank_work=exact_integer_rank(rows,context)
    facts={'dimension':rank,'rank_exact':rank,'ambient_dimension':d,'cardinality':n,'norm_squared':norm,
           'max_inner_product':maximum,'pairs_checked':n*(n-1)//2}
    return Verification('exact','integer coordinates with positive diagonal quadratic form and exact rank',
                        {**facts,'facts':facts,'data_sha256':digest(data),'integer_backend':'bounded exact binary32' if binary32 else ('int64 proven overflow-safe' if safe else 'Python integers'),
                         'rank_scalar_work':rank_work})


def compare_facts(claimed,computed):
    # Mathematical replay cannot certify novelty or unrestricted optimality.
    # Explicit negative flags are allowed on bounded search receipts.
    for key in ('new_sota','global_optimality_claim','worldwide_priority_claim'):
        if key in claimed and claimed[key] is not False:
            raise ValueError(f'{key} cannot be established by a construction verifier')
    protected = {'dimension','ambient_dimension','dimension_scope','rank_exact',
                 'rank_numerical','cardinality','length','min_distance','distance',
                 'norm_squared','max_inner_product','max_inner_product_upper_bound',
                 'max_cosine','max_cosine_upper_bound','finite_group_proved'}
    for key in set(claimed)&protected:
        if key not in computed:
            raise ValueError(f'Claimed {key} is not certified by this verifier')
    for key in set(claimed)&set(computed):
        if canonical_json(claimed[key])!=canonical_json(computed[key]):
            raise ValueError(f'Claimed {key} does not match verification: {claimed[key]} != {computed[key]}')


def verify_artifact(artifact,context=None):
    context=context or RunContext(Path('.'))
    if artifact.kind=='quadratic_spherical_code':
        from .engines.signed_layers import verify_artifact as verify_signed_layers
        report=verify_signed_layers(artifact,context)
    elif artifact.kind=='spherical_code' and artifact.data.get('backend')=='exact_integer_orbit':
        from .engines.exact_orbits import verify_artifact as verify_exact_orbit
        report=verify_exact_orbit(artifact,context)
    elif artifact.kind=='spherical_code' and artifact.data.get('backend')=='certified_coordinate_section':
        from .sections import verify_certified_section
        report=verify_certified_section(artifact,context)
    elif artifact.kind=='gram_spherical_code' and 'gluing' in artifact.data:
        from .compatibility import verify_glued_artifact
        report=verify_glued_artifact(artifact,context)
    elif artifact.kind=='gram_spherical_code':
        from .exact_gram import verify_rational_gram
        report=verify_rational_gram(artifact.data,context)
    elif artifact.kind=='spherical_code' and 'section' in artifact.data:
        from .engines.geometry import verify_section
        report=verify_section(artifact,context)
    elif artifact.kind=='syndrome_oracle' or (artifact.kind=='binary_code' and artifact.data.get('backend')=='binary_syndrome_cosets'):
        from .engines.cosets import verify_artifact as verify_coset
        report=verify_coset(artifact,context)
    elif artifact.kind=='spherical_code' and 'coordinates' in artifact.data:
        report=verify_integer_spherical(artifact.data,context)
    elif artifact.kind=='implicit_spherical_code' and artifact.data.get('backend')=='real_stabilizer_sector':
        from .engines.stabilizers import verify_artifact as verify_stabilizer
        report=verify_stabilizer(artifact,context)
    elif artifact.kind=='implicit_spherical_code' and artifact.data.get('backend')=='signed_support_layers':
        from .engines.layers import verify_artifact as verify_layers
        report=verify_layers(artifact,context)
    elif artifact.kind=='implicit_spherical_code':
        from .engines.geometry import verify_sign_embedding
        report=verify_sign_embedding(artifact,context)
    elif artifact.kind in ('binary_code','ring_code','polynomial_pair'):
        from .engines.rings import verify_artifact as verify_ring
        report=verify_ring(artifact,context)
    elif artifact.kind=='search_result' and artifact.data.get('backend')=='binary_translate_search':
        from .engines.binary_translate_search import verify_artifact as verify_binary_translate
        report=verify_binary_translate(artifact,context)
    elif artifact.kind=='search_result' and artifact.data.get('backend')=='gray_coset_assembly':
        from .engines.coset_search import verify_artifact as verify_coset_search
        report=verify_coset_search(artifact,context)
    elif artifact.kind in ('unitary_representation','spherical_code','search_result'):
        from .engines.groups import verify_artifact as verify_group
        report=verify_group(artifact,context)
    else:
        raise ValueError(f'No independent verifier for artifact kind {artifact.kind!r}')
    compare_facts(artifact.facts,report.details.get('facts',report.details))
    return report

