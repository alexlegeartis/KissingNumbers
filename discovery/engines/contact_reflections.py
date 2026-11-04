"""Verified engine adapter for the shared standalone contact-reflection code."""
import hashlib
from pathlib import Path

from ..core import Artifact, check_params, integer, digest
from ..certification import verify_artifact, verify_integer_spherical
from .construction_search import _seconds, blocker_exchange, integer_pool_graph
from .selection import maximum_independent_bipartite


def contact_reflections(inputs,params,context):
    from research import contact_reflection_20261104 as prototype
    check_params(params,('seeds','pair_limit','candidate_limit','blocker_limit','seconds'))
    if len(inputs)!=1 or inputs[0].kind!='spherical_code' or 'coordinates' not in inputs[0].data:
        raise ValueError('contact_reflections needs one explicit integer code')
    source=inputs[0];rows=source.data['coordinates']
    if not isinstance(rows,list) or not rows:raise ValueError('Nonempty source rows required')
    n=len(rows)
    cap=integer(params.get('candidate_limit',256),'candidate_limit')
    pairs=integer(params.get('pair_limit',100000),'pair_limit')
    seconds=_seconds(params.get('seconds',1),context)
    cutoff=params.get('blocker_limit')
    if cutoff is not None:integer(cutoff,'blocker_limit',0)
    seeds=params.get('seeds',list(range(min(n,200))))
    if not isinstance(seeds,list) or any(type(i) is not int or not 0<=i<n for i in seeds) or len(set(seeds))!=len(seeds):
        raise ValueError('Distinct valid integer seed indices required')
    context.limit('max_trials',len(seeds))
    context.limit('max_materialized_points',n+cap)
    context.limit('max_codewords',cap)
    context.limit('max_pairs',n*(n-1)//2+(n+cap)*(n+cap-1)//2+len(seeds)*n+n*cap+cap*(cap-1)//2)
    if verify_artifact(source,context).status!='exact':raise ValueError('Exact source required')
    d=len(rows[0]);norm=source.data['norm_squared'];weights=source.data.get('weights')
    if weights is None:weights=[1]*d
    context.limit('max_rank_work',n*d*d+(len(seeds)*n+n*cap+cap*cap)*d)
    metric=[[weights[i] if i==j else 0 for j in range(d)] for i in range(d)]
    pool=prototype.search(rows,metric,seeds=seeds,pair_limit=pairs,candidate_limit=cap,
                          blocker_limit=cutoff,seconds=seconds,return_pool=True)
    candidates,blocked=pool['candidates'],pool['blockers']
    receipt={**pool['receipt'],'seeds':seeds,'source_data_sha256':digest(source.data),
        'prototype_sha256':hashlib.sha256(Path(prototype.__file__).read_bytes()).hexdigest(),
        'source_validity':'independently verified in the engine',
        'metric_positive_definiteness':'positive diagonal metric independently checked',
        'new_sota':False}
    selected,removed=[],[]
    if candidates:
        matching=maximum_independent_bipartite(blocked,n)
        if len(matching['matching'])==len(candidates):
            receipt['matching_source_indices']=[j for _,j in matching['matching']]
            receipt['combinatorial_bound']='matching saturates this generated candidate pool'
        else:
            unique,_,left,right,original=integer_pool_graph(candidates,weights,norm,context)
            if len(unique)!=len(candidates):raise ValueError('Prototype returned duplicate candidates')
            proposal=blocker_exchange([],{'blockers':blocked,'candidate_conflicts':[list(e) for e in zip(left,right)],
                                          'seconds':seconds,'presolve':False},context)
            selected=proposal.data['selected_candidate_indices'];removed=proposal.data['removed_source_indices']
            receipt['solver']=proposal.data['solver']
    deleted=set(removed)
    receipt.update(gain=len(selected)-len(removed),selected_candidate_indices=selected,
                   removed_source_indices=removed,
                   scope='Finite contact/opposite-contact reflections and explicit generation/cost limits; exact source and complete output proof; no global optimality or priority claim')
    data={'coordinates':[row for i,row in enumerate(rows) if i not in deleted]+[candidates[i] for i in selected],
          'weights':weights,'norm_squared':norm,'search':receipt}
    report=verify_integer_spherical(data,context)
    report.details['contact_reflections']=receipt
    return Artifact('spherical_code',data,report.details['facts'],report)


OPERATORS={'search.contact_reflections':contact_reflections}
INPUT_KINDS={'search.contact_reflections':('spherical_code',)}
OUTPUT_KINDS={'search.contact_reflections':('spherical_code',)}
