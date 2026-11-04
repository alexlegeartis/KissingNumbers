"""Exact quadratic-average screening of full generated orthogonal orbits.

No conjugacy-class, irreducibility, finiteness or partial-orbit assumption is
silently certified. An explicit action set can be averaged and checked directly.
"""
from fractions import Fraction
from pathlib import Path

from .core import RunContext
from .exact_gram import bounded_fraction,check_fraction_size,rational_json


def analyze_orthogonal_actions(seed,actions,context=None):
    """Verify exact actions and screen a full orbit before enumerating it.

    Every listed action belongs to the group it generates, so one distinct
    violating image is an exact obstruction to the full generated orbit (and
    any larger orbit containing it). A selected partial orbit is not pruned.
    """
    context=context or RunContext(Path('.'))
    if not isinstance(seed,list) or not seed:
        raise ValueError('A nonzero rational seed is required')
    if not isinstance(actions,list) or not actions:
        raise ValueError('A nonempty exact action list is required')
    d,m=len(seed),len(actions)
    context.limit('max_dimension',d)
    context.limit('max_nodes',m)
    context.limit('max_gram_entries',m*d*d)
    context.limit('max_rank_work',3*m*d*d*d)
    if any(not isinstance(a,list) or len(a)!=d or any(not isinstance(row,list) or len(row)!=d for row in a) for a in actions):
        raise ValueError('Action matrices must match seed dimension')
    v=[bounded_fraction(x,context,'seed entry') for x in seed]
    def checked(value):return check_fraction_size(value,context)
    norm=checked(sum((checked(x*x) for x in v),Fraction(0)))
    if not norm:
        raise ValueError('A nonzero rational seed is required')
    total=[[Fraction(0) for _ in range(d)] for _ in range(d)]
    fixed=[];dots=[];witnesses=[]
    for index,raw in enumerate(actions):
        matrix=[[bounded_fraction(x,context,'action entry') for x in row] for row in raw]
        for i in range(d):
            for j in range(d):
                value=checked(sum((checked(matrix[k][i]*matrix[k][j]) for k in range(d)),Fraction(0)))
                if value!=int(i==j):
                    raise ValueError('Action is not exactly orthogonal')
                total[i][j]=checked(total[i][j]+checked(matrix[i][j]+matrix[j][i]))
        image=[checked(sum((checked(x*y) for x,y in zip(row,v)),Fraction(0))) for row in matrix]
        dot=checked(sum((checked(x*y) for x,y in zip(image,v)),Fraction(0)))
        dots.append(dot)
        if image==v:
            fixed.append(index)
        elif 2*dot>norm:
            witnesses.append({'action_index':index,'image':[rational_json(x) for x in image],
                              'dot':rational_json(dot),'norm_squared':rational_json(norm)})
    average=[[checked(x/(2*m)) for x in row] for row in total]
    scalar=average[0][0]
    scalar_average=all(average[i][j]==(scalar if i==j else 0) for i in range(d) for j in range(d))
    result={'dimension':d,'action_count':m,'norm_squared':rational_json(norm),
            'fixed_action_indices':fixed,'fixed_action_count':len(fixed),
            'normalized_dots':[rational_json(checked(x/norm)) for x in dots],
            'symmetric_average_is_scalar':scalar_average,
            'whole_generated_orbit_ruled_out':bool(witnesses),'violating_images':witnesses,
            'partial_orbits_pruned':False,'finite_group_proved':False,
            'conjugacy_class_proved':False,'irreducibility_proved':False,
            'scope':'Full orbit under the listed exact orthogonal actions; any partial-orbit selection remains outside this obstruction'}
    if scalar_average:
        required=checked(m*(2*scalar-1))
        minimum=max(0,-(-required.numerator//required.denominator))
        result.update(scalar_average=rational_json(scalar),
                      necessary_fixed_actions_for_kissing=minimum,
                      fixed_action_budget_violated=len(fixed)<minimum,
                      equality_constraints_apply=scalar==Fraction(1,2) and not fixed,
                      equality_constraints_satisfied=all(2*x==norm for x in dots))
        # A contradiction here would be a checker bug, not a mathematical result.
        if len(fixed)<minimum and not witnesses:
            raise ValueError('Average obstruction lacks its explicit violating image')
    return result
