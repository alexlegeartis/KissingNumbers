"""Joint shape assignment to character-code cosets.

An alternative is a whole coset of signed coordinate shapes. A Walsh transform
computes every cross product, including those inside the proposed families.
The optimizer may replace any incumbent family; its dual is never a proof.
"""
from time import monotonic

import numpy as np

from ..core import Artifact, check_params, integer
from ..certification import verify_integer_spherical
from .construction_search import _binary_solution, _diagnostic, _seconds


def _binary_rows(rows, width, name):
    if not isinstance(rows, list) or any(not isinstance(row, list) or len(row) != width
            or any(type(x) is not int or x not in (0, 1) for x in row) for row in rows):
        raise ValueError(f'{name} must be rectangular binary integer rows')
    pivots = {}
    for row in rows:
        value = sum(x << i for i, x in enumerate(row))
        while value:
            bit = value.bit_length()-1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
        if not value:
            raise ValueError(f'{name} must be linearly independent')
    return rows


def walsh_products(generators, shapes, weights, context):
    """Exact dot[a,b,m] = <shape[a], character(m)*shape[b]>.

    No full vector-pair sweep, assumed group action, or distance bound is used.
    Large coordinates use Python integers, including inside the transform.
    """
    if not isinstance(shapes, list) or not shapes or not isinstance(shapes[0], list) or not shapes[0]:
        raise ValueError('Nonempty coordinate shapes required')
    d, labels = len(shapes[0]), len(shapes)
    context.limit('max_dimension', d)
    context.limit('max_materialized_points', labels*d)
    if any(not isinstance(row, list) or len(row) != d or
           any(type(x) is not int for x in row) for row in shapes):
        raise ValueError('Shapes must be rectangular exact integer rows')
    if not isinstance(generators, list) or not generators:
        raise ValueError('Nonempty binary generators required')
    _binary_rows(generators, d, 'generators')
    size = 1 << len(generators)
    context.limit('max_codewords', size)
    context.limit('max_gram_entries', labels*labels*size)
    context.limit('max_rank_work', labels*labels*(size*len(generators)+d))
    if weights is None:
        weights = [1]*d
    if not isinstance(weights, list) or len(weights) != d:
        raise ValueError('Wrong character metric')
    for w in weights:
        integer(w, 'metric weight')
        context.limit('max_rational_bits', w.bit_length())
    for row in shapes:
        for x in row:
            context.limit('max_rational_bits', abs(x).bit_length())
    norms = [sum(w*x*x for w, x in zip(weights, row)) for row in shapes]
    if min(norms) <= 0 or len(set(norms)) != 1:
        raise ValueError('Shapes require a common positive squared norm')
    norm = norms[0]
    context.limit('max_rational_bits', norm.bit_length())
    columns = [sum(g[i] << j for j, g in enumerate(generators)) for i in range(d)]
    # Every intermediate butterfly is a signed sum of a subset of these
    # coefficients, so this absolute bound covers intermediates as well.
    bound = max(abs(x) for row in shapes for x in row)**2*sum(weights)
    dtype = np.int64 if bound < 2**62 else object
    dots = np.zeros((labels, labels, size), dtype=dtype)
    for a, first in enumerate(shapes):
        for b, second in enumerate(shapes):
            for column, w, x, y in zip(columns, weights, first, second):
                dots[a, b, column] += w*x*y
    step = 1
    while step < size:
        for start in range(0, size, 2*step):
            left = dots[:, :, start:start+step].copy()
            right = dots[:, :, start+step:start+2*step].copy()
            dots[:, :, start:start+step] = left+right
            dots[:, :, start+step:start+2*step] = left-right
        step *= 2
    return dots, weights, norm, columns


def family_problem(generators, shapes, weights, quotient, context):
    """Return exact forbidden-label tables and actual message cosets.

    One shape per coset is a search restriction. Every shape must have a free
    sign action; stabilizer quotients and unions of overlapping families are
    deliberately rejected rather than overcounted.
    """
    dots, weights, norm, columns = walsh_products(generators, shapes, weights, context)
    _binary_rows(quotient, len(generators), 'quotient')
    size, labels, qsize = dots.shape[2], len(shapes), 1 << len(quotient)
    context.limit('max_codewords', qsize*labels)
    masks = [sum(x << i for i, x in enumerate(row)) for row in quotient]
    syndromes = [sum(((m & mask).bit_count() & 1) << i for i, mask in enumerate(masks))
                 for m in range(size)]
    cosets = [[] for _ in range(qsize)]
    for m, q in enumerate(syndromes):
        cosets[q].append(m)
    forbidden = np.empty((qsize, labels, labels), dtype=bool)
    for q, messages in enumerate(cosets):
        forbidden[q] = 2*np.max(dots[:, :, messages], axis=2) > norm
    admissible = []
    for a in range(labels):
        if any(int(x) == norm for x in dots[a, a, 1:]):
            raise ValueError('Shape has a nontrivial sign stabilizer; free action required')
        admissible.append(not any(2*int(dots[a, a, m]) > norm for m in cosets[0] if m))
    return {'forbidden': forbidden, 'admissible': admissible, 'cosets': cosets,
            'weights': weights, 'norm': norm, 'columns': columns,
            'kernel_size': len(cosets[0]), 'message_count': size}


def check_assignment(problem, assignment):
    table, valid = problem['forbidden'], problem['admissible']
    n, labels = table.shape[:2]
    if not isinstance(assignment, list) or len(assignment) != n or any(
            type(a) is not int or not -1 <= a < labels for a in assignment):
        raise ValueError('Assignment needs one label or -1 per coset')
    for q, a in enumerate(assignment):
        if a < 0:
            continue
        if not valid[a]:
            raise ValueError('Selected character family has an internal conflict')
        for r in range(q):
            b = assignment[r]
            if b >= 0 and table[q ^ r, a, b]:
                raise ValueError('Character families have an exact cross conflict')
    return sum(a >= 0 for a in assignment)


def optimize_families(problem, context, seconds, initial=None, active=None):
    """Joint MILP, optionally a large neighborhood with a checked incumbent.

    Every active coset can change label or be deleted. No ordering of new
    families or acceptance of a solver tolerance enters the geometry proof.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    table, valid = problem['forbidden'], problem['admissible']
    n, labels = table.shape[:2]
    initial = [-1]*n if initial is None else initial[:]
    initial_count = check_assignment(problem, initial)
    if active is None:
        active = list(range(n))
    if not isinstance(active, list) or not active or len(set(active)) != len(active) or any(
            type(q) is not int or not 0 <= q < n for q in active):
        raise ValueError('active_cosets must be distinct coset indices')
    context.limit('max_codewords', len(active)*labels)
    active_set = set(active)
    frozen = [(q, a) for q, a in enumerate(initial) if q not in active_set and a >= 0]
    variables = [(q, a) for q in active for a in range(labels) if valid[a] and
                 not any(table[q ^ r, a, b] for r, b in frozen)]
    if not variables:
        return initial, {'status': 'no_admissible_variables', 'initial_families': initial_count}
    by_coset = {q: [] for q in active}
    for i, (q, a) in enumerate(variables):
        by_coset[q].append((i, a))
    rr, cc, row = [], [], 0
    for q in active:
        for i, a in by_coset[q]:
            rr.append(row); cc.append(i)
        row += 1
    # Cayley differences with no forbidden label pair need no graph entries.
    pairs = {d: np.argwhere(table[d]).tolist() for d in range(1, n) if table[d].any()}
    context.limit('max_rank_work', len(active)*sum(len(p) for p in pairs.values()))
    lookup = {(q, a): i for i, (q, a) in enumerate(variables)}
    edges = 0
    for q in active:
        for d, forbidden_pairs in pairs.items():
            r = q ^ d
            if r >= q or r not in active_set:
                continue
            for a, b in forbidden_pairs:
                i, j = lookup.get((q, a)), lookup.get((r, b))
                if i is not None and j is not None:
                    rr.extend((row, row)); cc.extend((i, j)); row += 1; edges += 1
                    context.limit('max_gram_entries', len(rr))
    matrix = coo_matrix((np.ones(len(rr)), (rr, cc)), shape=(row, len(variables))).tocsc()
    result = milp(-np.ones(len(variables)), integrality=np.ones(len(variables)),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, 0, 1),
                  options={'time_limit': seconds, 'mip_rel_gap': 0.0})
    proposal = _binary_solution(result, len(variables))
    best = initial
    if proposal is not None:
        candidate = initial[:]
        for q in active:
            candidate[q] = -1
        for i in np.flatnonzero(proposal):
            q, a = variables[int(i)]
            if candidate[q] >= 0:
                raise ValueError('Solver selected two shapes in the same coset')
            candidate[q] = a
        if check_assignment(problem, candidate) > initial_count:
            best = candidate
    diagnostic = {**_diagnostic(result), 'initial_families': initial_count,
                  'selected_families': check_assignment(problem, best),
                  'active_cosets': len(active), 'variables': len(variables), 'conflicts': edges}
    return best, diagnostic


def materialize_assignment(problem, shapes, assignment):
    check_assignment(problem, assignment)
    rows = []
    for q, a in enumerate(assignment):
        if a >= 0:
            for m in problem['cosets'][q]:
                rows.append([x*(1-2*((m & col).bit_count() & 1))
                             for x, col in zip(shapes[a], problem['columns'])])
    return rows


def search(inputs, params, context):
    check_params(params, ('generators', 'shapes', 'weights', 'quotient', 'seconds',
                          'initial_assignment', 'active_cosets'))
    if inputs:
        raise ValueError('character_families takes no inputs')
    seconds = _seconds(params.get('seconds', 1), context)
    started = monotonic()
    problem = family_problem(params.get('generators'), params.get('shapes'), params.get('weights'),
                             params.get('quotient', []), context)
    context.limit('max_materialized_points', problem['message_count'])
    context.limit('max_pairs', problem['message_count']*(problem['message_count']-1)//2)
    assignment, diagnostic = optimize_families(problem, context, seconds,
        params.get('initial_assignment'), params.get('active_cosets'))
    rows = materialize_assignment(problem, params['shapes'], assignment)
    if not rows:
        raise ValueError('No nonempty feasible character-family construction found')
    data = {'coordinates': rows, 'weights': problem['weights'], 'norm_squared': problem['norm'],
            'search': {'method': 'Walsh character kernel and joint coset-shape assignment',
                       'assignment': assignment, 'kernel_size': problem['kernel_size'],
                       'solver': diagnostic, 'seconds': monotonic()-started,
                       'target_counts_used': False,
                       'scope': 'One shape per coset, free sign action; numerical optimization is not an upper bound'}}
    report = verify_integer_spherical(data, context)
    return Artifact('spherical_code', data, report.details['facts'], report)


def assemble_orbits(inputs, params, context):
    """Select arbitrary disjoint sign orbits, including sparse stabilizers.

    This variant allows many components in the same sign-code coset. Thus a
    root/octad base and several cap families can all be deleted and replaced
    together. Distinct overlapping orbits conflict rather than double count.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    check_params(params, ('generators', 'shapes', 'weights', 'seconds', 'initial_shapes'))
    if inputs:
        raise ValueError('character_orbits takes no inputs')
    seconds = _seconds(params.get('seconds', 1), context)
    shapes = params.get('shapes')
    dots, weights, norm, columns = walsh_products(params.get('generators'), shapes,
                                                 params.get('weights'), context)
    labels, _, size = dots.shape
    valid, counts = [], []
    for a in range(labels):
        stabilizers = sum(int(x) == norm for x in dots[a, a])
        if not stabilizers or size % stabilizers:
            raise ValueError('Invalid sign stabilizer count')
        counts.append(size//stabilizers)
        valid.append(not any(2*int(x) > norm and int(x) != norm for x in dots[a, a]))
    bad = 2*np.max(dots, axis=2) > norm
    initial = params.get('initial_shapes', [])
    if not isinstance(initial, list) or len(set(initial)) != len(initial) or any(
            type(a) is not int or not 0 <= a < labels or not valid[a] for a in initial):
        raise ValueError('initial_shapes must be distinct internally compatible orbit indices')
    if any(bad[a,b] for i,a in enumerate(initial) for b in initial[:i]):
        raise ValueError('Initial character orbits conflict or overlap')
    variables = [a for a in range(labels) if valid[a]]
    potential = sum(counts[a] for a in variables)
    # Allocation bound uses the union of all admissible proposed families,
    # before optimization. Users can explicitly increase it for large pools.
    context.limit('max_materialized_points', potential)
    context.limit('max_pairs', potential*(potential-1)//2)
    edges = [(i,j) for i,a in enumerate(variables) for j,b in enumerate(variables[:i]) if bad[a,b]]
    context.limit('max_gram_entries', 2*len(edges))
    rr,cc=[],[]
    for r,(i,j) in enumerate(edges):
        rr.extend((r,r));cc.extend((i,j))
    best, diagnostic = initial[:], {'status':'no_admissible_orbits'}
    if variables:
        matrix = coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(edges),len(variables))).tocsc()
        result = milp(-np.array([counts[a] for a in variables],float),
                      integrality=np.ones(len(variables)), bounds=Bounds(0,1),
                      constraints=LinearConstraint(matrix,0,1),
                      options={'time_limit':seconds,'mip_rel_gap':0.0})
        proposal = _binary_solution(result,len(variables))
        if proposal is not None:
            candidate=[variables[int(i)] for i in np.flatnonzero(proposal)]
            if any(bad[a,b] for i,a in enumerate(candidate) for b in candidate[:i]):
                raise ValueError('Rounded character orbit solution conflicts')
            if sum(counts[a] for a in candidate)>sum(counts[a] for a in best):
                best=candidate
        diagnostic=_diagnostic(result)
    rows=[]
    for a in best:
        orbit={tuple(x*(1-2*((m&col).bit_count()&1)) for x,col in zip(shapes[a],columns))
               for m in range(size)}
        if len(orbit)!=counts[a]:
            raise ValueError('Stabilizer orbit count does not match materialization')
        rows.extend(map(list,sorted(orbit)))
    if not rows:
        raise ValueError('No nonempty character orbit assembly found')
    data={'coordinates':rows,'weights':weights,'norm_squared':norm,
          'search':{'method':'Joint weighted character-orbit assembly with exact stabilizers',
                    'selected_shapes':best,'initial_points':sum(counts[a] for a in initial),
                    'admissible_orbits':len(variables),'proposed_orbits':labels,
                    'conflicts':len(edges),'solver':diagnostic,'target_counts_used':False,
                    'all_initial_families_removable':True,
                    'scope':'Finite proposed whole sign orbits; no optimality or novelty proof'}}
    report=verify_integer_spherical(data,context)
    return Artifact('spherical_code',data,report.details['facts'],report)


OPERATORS = {'search.character_families': search, 'search.character_orbits': assemble_orbits}
INPUT_KINDS = {op: () for op in OPERATORS}
OUTPUT_KINDS = {op: ('spherical_code',) for op in OPERATORS}
