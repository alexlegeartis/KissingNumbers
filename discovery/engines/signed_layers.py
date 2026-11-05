"""Conditional signed-code switching and exact Q(sqrt(3)) coupled layers.

All output points are independently checked; no antipodality is required.
Positive integer squared norms are normalized per point. Raw coordinates lie
in Q(sqrt(3)); normalized coordinates may require further square roots.
"""
from pathlib import Path
from collections import Counter
import numpy as np
from ..core import Artifact, Verification, RunContext, check_params, integer, digest


def require(condition, message):
    if not condition: raise ValueError(message)


def sign_sqrt3(a, b):
    """Exact sign of a+b*sqrt(3), for Python integers."""
    if not b: return (a > 0)-(a < 0)
    if not a: return (b > 0)-(b < 0)
    if a > 0 and b > 0: return 1
    if a < 0 and b < 0: return -1
    difference = a*a-3*b*b
    return ((difference > 0)-(difference < 0)) * ((a > 0)-(a < 0))


def compatible_dot(a, b, norm_product):
    """Normalize two raw vectors: dot <= sqrt(norm_product)/2."""
    if sign_sqrt3(a, b) <= 0: return True
    return sign_sqrt3(4*(a*a+3*b*b)-norm_product, 8*a*b) <= 0


def _full_rank_witness(A, B, context):
    # sqrt(3)->4 modulo13 is a ring homomorphism (4^2=3 mod13).
    # A nonzero minor modulo13 proves that exact algebraic minor nonzero.
    n, d = len(A), len(A[0]);context.limit('max_rank_work', n*d*d)
    pivots = {}
    for a, b in zip(A, B):
        row = [(x+4*y) % 13 for x, y in zip(a,b)]
        for pivot, basis in pivots.items():
            factor = row[pivot]
            if factor: row = [(x-factor*y) % 13 for x,y in zip(row,basis)]
        index = next((i for i,x in enumerate(row) if x), None)
        if index is not None:
            inverse = pow(row[index], -1, 13)
            pivots[index] = [(x*inverse) % 13 for x in row]
            if len(pivots) == d: return d
    return len(pivots)


def verify_artifact(artifact, context=None):
    context = context or RunContext(Path('.'))
    require(artifact.kind == 'quadratic_spherical_code', 'Expected quadratic_spherical_code')
    data = artifact.data
    require(set(data) == {'radicand','rational_coordinates','radical_coordinates'} and data['radicand'] == 3 and type(data['radicand']) is int, 'Unsupported quadratic coordinate fields')
    A, B = data['rational_coordinates'], data['radical_coordinates']
    require(isinstance(A,list) and A and isinstance(A[0],list) and A[0], 'Nonempty coordinate matrix required')
    n,d = len(A),len(A[0]);context.limit('max_materialized_points',n);context.limit('max_dimension',d)
    context.limit('max_pairs',n*(n-1)//2);context.limit('max_rank_work',n*d*d)
    require(isinstance(B,list) and len(B)==n, 'Radical coordinate row count')
    for row in A+B:
        require(isinstance(row,list) and len(row)==d, 'Ragged quadratic coordinates')
        require(all(type(x) is int for x in row), 'Quadratic coordinate coefficients must be integers')
        for x in row: context.limit('max_rational_bits',abs(x).bit_length())
    norms=[]
    for a,b in zip(A,B):
        require(sum(x*y for x,y in zip(a,b))==0, 'Only rational squared norms are supported')
        norm=sum(x*x+3*y*y for x,y in zip(a,b));require(norm>0,'Zero point')
        norms.append(norm)
    maximum = max(abs(x) for row in A+B for x in row)
    dtype = np.int64 if 4*d*maximum*maximum < 2**62 else object
    AA,BB=np.array(A,dtype=dtype),np.array(B,dtype=dtype)
    # When both algebraic embeddings have the same rational norm, exact
    # Cauchy bounds keep every squared sign comparison below 243*max(norm)^4.
    # Larger inputs retain the arbitrary-precision scalar fallback.
    vectorized=dtype is np.int64 and 243*max(norms)**4<2**62
    if vectorized:
        NN=np.asarray(norms,dtype=np.int64)
        def positive(p,q):
            return ((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q))
    for start in range(0,n,96):
        P=AA[start:start+96]@AA.T+3*(BB[start:start+96]@BB.T)
        Q=AA[start:start+96]@BB.T+BB[start:start+96]@AA.T
        if vectorized:
            R=4*(P*P+3*Q*Q)-NN[start:start+len(P),None]*NN[None,:];S=8*P*Q
            bad=positive(P,Q)&positive(R,S)
            for offset in range(len(P)):bad[offset,:start+offset+1]=False
            require(not np.any(bad), 'Quadratic normalized pair violation or duplicate')
        else:
            for offset in range(len(P)):
                i=start+offset
                for j in range(i+1,n):
                    require(compatible_dot(int(P[offset,j]),int(Q[offset,j]),norms[i]*norms[j]), 'Quadratic normalized pair violation or duplicate')
    rank_lower=_full_rank_witness(A,B,context)
    facts={'dimension':d,'ambient_dimension':d,'cardinality':n,
           'max_cosine_upper_bound':[1,2], 'rank_lower_bound':rank_lower,
           'dimension_scope':'ambient; modular minor lower bound'}
    if rank_lower==d:
        facts.update(rank_exact=d,dimension_scope='exact; nonzero algebraic minor proved modulo13')
    if len(set(norms))==1: facts['norm_squared']=norms[0]
    from ..certification import compare_facts
    compare_facts(artifact.facts,facts)
    return Verification('exact','all normalized pair inequalities in Q(sqrt(3)); positive rational squared norms',
                        {'facts':facts,'pairs_checked':n*(n-1)//2,'norm_histogram':{str(k):v for k,v in Counter(norms).items()},
                         'rank_lower_bound':rank_lower,'data_sha256':digest(data),
                         'arithmetic':'integer sign comparisons and exact squaring in Q(sqrt(3))',
                         'antipodality_required':False})


def conditional_swap(inputs,params,context):
    """Apply an explicit signed-coordinate predicate, then independently replay."""
    from ..certification import verify_artifact as verify, verify_integer_spherical
    check_params(params,('swap','conditions','support_masks'))
    require(len(inputs)==1 and inputs[0].kind=='spherical_code','An explicit spherical code is required')
    context=context or RunContext(Path('.'))
    source=inputs[0];require('coordinates' in source.data,'Explicit coordinates required')
    require(verify(source,context).status=='exact','Exact source required')
    rows=source.data['coordinates'];d=len(rows[0]);swap=params.get('swap')
    require(isinstance(swap,list) and len(swap)==2 and all(type(x) is int and 0<=x<d for x in swap) and swap[0]!=swap[1],'Two distinct valid swap indices required')
    weights=source.data.get('weights',[1]*d)
    require(weights[swap[0]]==weights[swap[1]],'Swap must preserve the supplied metric')
    conditions=params.get('conditions');require(isinstance(conditions,list) and conditions,'Conditions required')
    for condition in conditions:
        require(isinstance(condition,dict) and condition,'Condition must be a nonempty index/value object')
        require(all(isinstance(k,str) and k.isdigit() and 0<=int(k)<d and type(v) is int for k,v in condition.items()),'Invalid signed predicate')
    masks=params.get('support_masks')
    if masks is not None:require(isinstance(masks,list) and all(type(x) is int and 0<x<1<<d for x in masks),'Invalid support masks')
    output=[]
    for row in rows:
        changed=list(row);mask=sum(1<<i for i,x in enumerate(row) if x)
        if (masks is None or mask in masks) and any(all(row[int(k)]==v for k,v in condition.items()) for condition in conditions):
            changed[swap[0]],changed[swap[1]]=changed[swap[1]],changed[swap[0]]
        output.append(changed)
    data={'coordinates':output,'weights':weights,'norm_squared':source.data['norm_squared']}
    report=verify_integer_spherical(data,context)
    return Artifact('spherical_code',data,report.details['facts'],report)


def coupled_lift(inputs,params,context):
    """Lift an exact integer equator and supplied Q(sqrt3) cap bodies to layers.

    Cap bodies, heights, pole amplitude and equator scale are parameters. Every
    cross pair is checked; no compatibility is inferred merely from the recipe.
    """
    from ..certification import verify_artifact as verify
    check_params(params,('equator_scale','cap_rational','cap_radical','heights','pole_amplitude'))
    require(len(inputs)==1 and inputs[0].kind=='spherical_code','One explicit integer equator is required')
    context=context or RunContext(Path('.'))
    source=inputs[0];require('coordinates' in source.data,'Explicit equator required')
    require(verify(source,context).status=='exact','Exact source required')
    rows=source.data['coordinates'];d=len(rows[0]);require(source.data.get('weights',[1]*d)==[1]*d,'Euclidean equator required')
    scale=integer(params.get('equator_scale',1),'equator_scale')
    context.limit('max_rational_bits',scale.bit_length())
    caps=params.get('cap_rational');rad=params.get('cap_radical')
    require(isinstance(caps,list) and isinstance(rad,list) and len(caps)==len(rad),'Cap body matrices required')
    heights=params.get('heights');require(isinstance(heights,list) and heights,'Integer heights required')
    total=len(rows)+len(caps)*len(heights)+2
    context.limit('max_materialized_points',total);context.limit('max_dimension',d+1)
    context.limit('max_pairs',total*(total-1)//2);context.limit('max_rank_work',total*(d+1)*(d+1))
    require(all(isinstance(row,list) and len(row)==d and all(type(x) is int for x in row) for row in caps+rad),'Invalid cap body coefficients')
    require(all(type(h) is int for h in heights) and len(set(heights))==len(heights),'Distinct integer heights required')
    pole=integer(params.get('pole_amplitude'),'pole_amplitude')
    context.limit('max_rational_bits',pole.bit_length())
    for row in caps+rad:
        for x in row:context.limit('max_rational_bits',abs(x).bit_length())
    for h in heights:context.limit('max_rational_bits',abs(h).bit_length())
    source_max=max(abs(x) for row in rows for x in row)
    context.limit('max_rational_bits',(source_max*scale).bit_length())
    A=[[scale*x for x in row]+[0] for row in rows];B=[[0]*(d+1) for row in rows]
    for height in heights:
        A.extend([row+[height] for row in caps]);B.extend([row+[0] for row in rad])
    A.extend([[0]*d+[pole],[0]*d+[-pole]]);B.extend([[0]*(d+1),[0]*(d+1)])
    artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':A,'radical_coordinates':B})
    artifact.verification=verify_artifact(artifact,context);artifact.facts=artifact.verification.details['facts']
    return artifact


def free_zone_join(inputs,params,context):
    """Join two signed weight-four codes with caller-supplied coloured edges.

    Two-coordinate edges from equally indexed colour classes are paired and
    lifted to all sixteen signs. Every resultant pair is independently checked.
    """
    from itertools import product
    from ..certification import verify_artifact as verify, verify_integer_spherical
    check_params(params,('left_matchings','right_matchings'))
    context=context or RunContext(Path('.'))
    require(len(inputs)==2 and all(a.kind=='spherical_code' and 'coordinates' in a.data for a in inputs),'Two explicit signed codes required')
    for source in inputs:
        require(verify(source,context).status=='exact','Exact signed-code inputs required')
        rows=source.data['coordinates'];d=len(rows[0])
        require(source.data.get('weights',[1]*d)==[1]*d and source.data['norm_squared']==4,'Euclidean norm-four input codes required')
        require(all(all(x in (-1,0,1) for x in row) for row in rows),'Signed incidence coordinates required')
    da,db=(len(a.data['coordinates'][0]) for a in inputs)
    left,right=params.get('left_matchings'),params.get('right_matchings')
    require(isinstance(left,list) and isinstance(right,list) and len(left)==len(right) and left,'Paired colour classes required')
    for groups,d in ((left,da),(right,db)):
        for group in groups:
            require(isinstance(group,list) and group,'Nonempty edge colour required')
            for edge in group:
                require(isinstance(edge,list) and len(edge)==2 and all(type(i) is int and 0<=i<d for i in edge) and edge[0]!=edge[1],'Invalid coloured edge')
    total=len(inputs[0].data['coordinates'])+len(inputs[1].data['coordinates'])+16*sum(len(a)*len(b) for a,b in zip(left,right))
    context.limit('max_dimension',da+db);context.limit('max_materialized_points',total);context.limit('max_pairs',total*(total-1)//2)
    output=[row+[0]*db for row in inputs[0].data['coordinates']]+[[0]*da+row for row in inputs[1].data['coordinates']]
    for first,second in zip(left,right):
        for a in first:
            for b in second:
                support=a+[i+da for i in b]
                for signs in product((-1,1),repeat=4):
                    row=[0]*(da+db)
                    for i,sign in zip(support,signs):row[i]=sign
                    output.append(row)
    data={'coordinates':output,'norm_squared':4}
    report=verify_integer_spherical(data,context)
    return Artifact('spherical_code',data,report.details['facts'],report)


OPERATORS={'signed_layers.conditional_swap':conditional_swap,'signed_layers.coupled_lift':coupled_lift,'signed_layers.free_zone_join':free_zone_join}
INPUT_KINDS={'signed_layers.conditional_swap':('spherical_code',),'signed_layers.coupled_lift':('spherical_code',),'signed_layers.free_zone_join':('spherical_code','spherical_code')}
