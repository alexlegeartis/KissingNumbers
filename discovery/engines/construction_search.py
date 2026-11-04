"""Target-blind support synthesis and exact finite point-pool assembly.

Numerical MILP is a proposal mechanism. Only the independently checked output
coordinates establish a lower bound; solver duals never establish exclusion.
"""
from itertools import combinations, product
from math import comb, isfinite
from time import monotonic

import numpy as np

from ..core import Artifact, check_params, digest, integer
from ..certification import verify_artifact, verify_integer_spherical


def _seconds(value, context):
    if type(value) not in (int, float) or not isfinite(value) or value <= 0:
        raise ValueError('seconds must be positive and finite')
    context.limit('solver_seconds', value)
    return float(value)


def _diagnostic(result):
    def finite(value):
        return float(value) if value is not None and isfinite(float(value)) else None
    return {'status': int(result.status), 'message': str(result.message),
            'numerical_dual_bound': finite(getattr(result, 'mip_dual_bound', None)),
            'numerical_gap': finite(getattr(result, 'mip_gap', None)),
            'scope': 'Numerical diagnostics only; no exact finite or global upper bound'}


def _binary_solution(result, count):
    if result.x is None:
        return None
    if len(result.x) != count or not np.all(np.isfinite(result.x)):
        raise ValueError('Malformed solver incumbent')
    values = np.rint(result.x).astype(np.int64)
    if np.any((values != 0) & (values != 1)) or np.max(abs(values-result.x), initial=0) > 1e-5:
        raise ValueError('Solver incumbent is not binary')
    return values


def _permutation(cycles, dimension):
    if not isinstance(cycles, list) or not cycles:
        raise ValueError('cycles must be a nonempty list of positive lengths')
    for value in cycles:
        integer(value, 'cycle length')
    if sum(cycles) != dimension:
        raise ValueError('Cycle lengths must sum to dimension')
    permutation, offset = [], 0
    for size in cycles:
        permutation.extend(list(range(offset+1, offset+size))+[offset])
        offset += size
    return permutation


def support_packing(inputs, params, context):
    """Choose weight-four supports with each triple used at most once.

    Optional permutation symmetry compresses variables into cyclic orbits.
    All 16 sign choices and the 2*n coordinate axes are constructed afresh.
    There is no table of target counts or supplied support designs.
    """
    from scipy.optimize import milp, Bounds, LinearConstraint
    from scipy.sparse import coo_matrix
    check_params(params, ('dimension', 'cycles', 'seconds'))
    if inputs:
        raise ValueError('support_packing takes no inputs')
    n = integer(params.get('dimension'), 'dimension')
    context.limit('max_dimension', n)
    seconds = _seconds(params.get('seconds', 1), context)
    count = comb(n, 4) if n >= 4 else 0
    context.limit('max_codewords', count)
    # Any packing uses each of C(n,3) triples at most once. This generic
    # counting bound preflights output allocation, independent of the solver.
    upper = comb(n, 3)//4 if n >= 4 else 0
    context.limit('max_materialized_points', 16*upper+2*n)
    context.limit('max_pairs', (16*upper+2*n)*(16*upper+2*n-1)//2)
    context.limit('max_gram_entries', 4*count)
    permutation = _permutation(params.get('cycles', [1]*n), n)
    # Orbit closure visits every support once; bit/operation count is bounded.
    context.limit('max_rank_work', max(1, count)*max(n, 4)*4)
    remaining = set(combinations(range(n), 4))
    orbits, invoices = [], []
    while remaining:
        current, orbit = min(remaining), set()
        while current not in orbit:
            orbit.add(current)
            current = tuple(sorted(permutation[i] for i in current))
        remaining.difference_update(orbit)
        triples = [t for block in orbit for t in combinations(block, 3)]
        if len(set(triples)) == len(triples):
            orbits.append(sorted(orbit))
            invoices.append(triples)
    started = monotonic()
    # A deterministic greedy packing supplies a feasible fallback even when a
    # solver times out before producing an incumbent.
    selected, used = [], set()
    for j in sorted(range(len(orbits)), key=lambda j: (-len(orbits[j]), orbits[j])):
        if not any(t in used for t in invoices[j]):
            selected.append(j)
            used.update(invoices[j])
    diagnostic = {'status': 'empty_support_domain'}
    if orbits:
        triples = {t:i for i,t in enumerate(combinations(range(n), 3))}
        rr, cc = [], []
        for j, invoice in enumerate(invoices):
            rr.extend(triples[t] for t in invoice)
            cc.extend([j]*len(invoice))
        matrix = coo_matrix((np.ones(len(rr)), (rr, cc)), shape=(len(triples), len(orbits))).tocsc()
        result = milp(-np.array(list(map(len, orbits)), dtype=float), integrality=np.ones(len(orbits)),
                      bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, 0, 1),
                      options={'time_limit': seconds, 'mip_rel_gap': 0.0})
        diagnostic = _diagnostic(result)
        proposal = _binary_solution(result, len(orbits))
        if proposal is not None:
            indices = np.flatnonzero(proposal).tolist()
            # Discrete incidence is checked with Python sets, not an LP tolerance.
            seen = set()
            for j in indices:
                if any(t in seen for t in invoices[j]):
                    raise ValueError('Rounded packing has a repeated triple')
                seen.update(invoices[j])
            if sum(len(orbits[j]) for j in indices) > sum(len(orbits[j]) for j in selected):
                selected = indices
    supports = sorted(block for j in selected for block in orbits[j])
    rows = []
    for i in range(n):
        for sign in (-2, 2):
            row = [0]*n
            row[i] = sign
            rows.append(row)
    for support in supports:
        for signs in product((-1, 1), repeat=4):
            row = [0]*n
            for i, sign in zip(support, signs):
                row[i] = sign
            rows.append(row)
    data = {'coordinates': rows, 'weights': [1]*n, 'norm_squared': 4,
            'search': {'method': 'triple-disjoint support orbit packing',
                       'cycles': params.get('cycles', [1]*n), 'supports': [list(s) for s in supports],
                       'support_domain': count, 'admissible_orbits': len(orbits),
                       'solver': diagnostic, 'seconds': monotonic()-started,
                       'target_counts_used': False}}
    report = verify_integer_spherical(data, context)
    return Artifact('spherical_code', data, report.details['facts'], report)


def integer_pool_graph(rows, weights, norm, context):
    """Deduplicate equal-norm rows; check every conflict exactly in blocks."""
    integer(norm, 'norm_squared')
    if not isinstance(rows, list) or not rows or not isinstance(rows[0], list) or not rows[0]:
        raise ValueError('Nonempty pool coordinates required')
    n, d = len(rows), len(rows[0])
    context.limit('max_materialized_points', n)
    context.limit('max_dimension', d)
    context.limit('max_pairs', n*(n-1)//2)
    context.limit('max_gram_entries', min(64, n)*n)
    if weights is None:
        weights = [1]*d
    if not isinstance(weights, list) or len(weights) != d:
        raise ValueError('Wrong pool metric')
    for w in weights:
        integer(w, 'pool metric weight')
        context.limit('max_rational_bits', w.bit_length())
    context.limit('max_rational_bits', norm.bit_length())
    for row in rows:
        if not isinstance(row, list) or len(row) != d or any(type(x) is not int for x in row):
            raise ValueError('Pool rows must be rectangular exact integers')
        for value in row:
            context.limit('max_rational_bits', abs(value).bit_length())
        if sum(w*x*x for w,x in zip(weights, row)) != norm:
            raise ValueError('Incorrect candidate norm')
    unique, lookup, original_indices = [], {}, []
    for row in rows:
        key = tuple(row)
        if key not in lookup:
            lookup[key] = len(unique)
            unique.append(row[:])
        original_indices.append(lookup[key])
    bound = max(abs(x) for row in unique for x in row)**2*sum(weights)
    safe = bound < 2**62 and max(weights) < 2**62
    dtype = np.int64 if safe else object
    values, metric = np.array(unique, dtype=dtype), np.array(weights, dtype=dtype)
    left, right = [], []
    for at in range(0, len(unique), 64):
        dots = (values[at:at+64]*metric)@values.T
        for offset, row in enumerate(dots):
            i = at+offset
            js = np.flatnonzero(row[i+1:] > norm//2)+i+1
            context.limit('max_gram_entries', 2*(len(left)+len(js)))
            left.extend([i]*len(js))
            right.extend(map(int, js))
    return unique, weights, left, right, original_indices


def select_integer_pool(rows, weights, norm, context, *, initial_indices=(), seconds=1):
    """Point-level selection with all initial points freely removable.

    Arbitrarily many mutually conflicting cohorts may be present. No whole
    orbit, antipodal pair, fiber, component or fixed core is enforced.
    """
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix
    seconds = _seconds(seconds, context)
    initial_indices = list(initial_indices)
    if any(type(i) is not int or not 0 <= i < len(rows) for i in initial_indices):
        raise ValueError('Initial pool indices out of range')
    unique, metric, left, right, indices = integer_pool_graph(rows, weights, norm, context)
    initial = set(indices[i] for i in initial_indices)
    if any(i in initial and j in initial for i,j in zip(left, right)):
        raise ValueError('Initial subset has an exact conflict')
    n = len(unique)
    neighbors = [set() for _ in unique]
    for i,j in zip(left, right):
        neighbors[i].add(j)
        neighbors[j].add(i)
    fallback = set(initial)
    for i in sorted(range(n), key=lambda i: (len(neighbors[i]), i)):
        if not (neighbors[i] & fallback):
            fallback.add(i)
    if left:
        # z_i deletes an initial point and inserts any other point. This is a
        # change of variables only; every point remains free to move.
        signs = np.array([-1 if i in initial else 1 for i in range(n)])
        cols = np.column_stack((left, right)).ravel()
        rr = np.repeat(np.arange(len(left)), 2)
        matrix = coo_matrix((signs[cols].astype(float), (rr, cols)), shape=(len(left), n)).tocsc()
        upper = np.array([1-int(i in initial)-int(j in initial) for i,j in zip(left,right)], dtype=float)
        result = milp(-signs.astype(float), integrality=np.ones(n), bounds=Bounds(0,1),
                      constraints=LinearConstraint(matrix, -np.inf, upper),
                      options={'time_limit': seconds, 'mip_rel_gap': 0.0})
        diagnostic = _diagnostic(result)
        proposed = _binary_solution(result, n)
        if proposed is not None:
            selected = {i for i in range(n) if (i in initial) != bool(proposed[i])}
            if any(i in selected and j in selected for i,j in zip(left,right)):
                raise ValueError('Rounded pool incumbent has an exact conflict')
            if len(selected) > len(fallback):
                fallback = selected
    else:
        diagnostic = {'status': 'edgeless_pool'}
    selected = sorted(fallback)
    data = {'coordinates': [unique[i] for i in selected], 'weights': metric, 'norm_squared': norm}
    report = verify_integer_spherical(data, context)
    optimization = {'pool_points_before_dedup': len(rows), 'distinct_pool_points': n,
                    'conflict_edges': len(left), 'pairs_checked': n*(n-1)//2,
                    'pool_sha256': digest({'coordinates': unique, 'weights': metric, 'norm_squared': norm}),
                    'selected_indices': selected, 'initial_indices': sorted(initial),
                    'all_initial_points_removable': True, 'solver': diagnostic,
                    'scope': 'Exact output geometry; numerical selection has no exact optimality claim'}
    data['search'] = optimization
    # The extra diagnostic data does not alter any mathematical fact.
    report = verify_integer_spherical(data, context)
    return Artifact('spherical_code', data, report.details['facts'], report)


def integer_pool(inputs, params, context):
    check_params(params, ('coordinates', 'weights', 'norm_squared', 'initial_indices', 'seconds'))
    if inputs:
        raise ValueError('integer_pool takes candidate coordinates as parameters')
    return select_integer_pool(params.get('coordinates'), params.get('weights'), params.get('norm_squared'), context,
                               initial_indices=params.get('initial_indices', ()), seconds=params.get('seconds', 1))


def contact_exchange(inputs, params, context):
    """Generate x+y at exact contact <x,y>=-norm/2, then assemble all points."""
    check_params(params, ('candidate_limit', 'seconds'))
    if len(inputs) != 1 or inputs[0].kind != 'spherical_code' or 'coordinates' not in inputs[0].data:
        raise ValueError('contact_exchange needs one explicit integer code')
    source = inputs[0]
    rows = source.data['coordinates']
    if not isinstance(rows, list) or not rows:
        raise ValueError('Nonempty source required')
    n, d = len(rows), len(rows[0])
    limit = integer(params.get('candidate_limit', 1000), 'candidate_limit', minimum=0)
    # Reject the possible full union before source replay or generation.
    context.limit('max_materialized_points', n+limit)
    context.limit('max_pairs', (n+limit)*(n+limit-1)//2)
    context.limit('max_dimension', d)
    if verify_artifact(source, context).status != 'exact':
        raise ValueError('Exact source required')
    norm = source.data['norm_squared']
    weights = source.data.get('weights', [1]*d)
    known = set(map(tuple, rows))
    candidates, origins, complete = [], [], True
    for i, a in enumerate(rows):
        for j in range(i+1, n):
            b = rows[j]
            if 2*sum(w*x*y for w,x,y in zip(weights,a,b)) != -norm:
                continue
            point = tuple(x+y for x,y in zip(a,b))
            if point in known:
                continue
            if len(candidates) == limit:
                complete = False
                break
            known.add(point)
            candidates.append(list(point))
            origins.append([i,j])
        if not complete:
            break
    result = select_integer_pool(rows+candidates, weights, norm, context,
                                 initial_indices=range(n), seconds=params.get('seconds', 1))
    result.data['search'].update(generation='exact negative-contact sums x+y', generated_candidates=len(candidates),
                                 generation_complete=complete, candidate_origins=origins)
    result.verification = verify_integer_spherical(result.data, context)
    return result


OPERATORS = {'search.support_packing': support_packing, 'selection.integer_pool': integer_pool,
             'search.contact_exchange': contact_exchange}
INPUT_KINDS = {'search.support_packing': (), 'selection.integer_pool': (),
               'search.contact_exchange': ('spherical_code',)}
OUTPUT_KINDS = {name: ('spherical_code',) for name in OPERATORS}
