"""Individual-point exchanges with a diagonal sign-group compatibility oracle.

Source symmetry accelerates screening; it is never imposed on the output.
Stabilizers are deduplicated, and internally conflicting orbits remain eligible.
"""
import numpy as np

from ..core import Artifact, check_params, digest, integer
from ..certification import verify_artifact, verify_integer_spherical
from .characters import _binary_rows
from .construction_search import blocker_exchange, _seconds


def sign_group(generators, dimension, context):
    """Message order respects multiplication by XOR, including the identity."""
    if not isinstance(generators, list):
        raise ValueError('generators must be binary rows')
    _binary_rows(generators, dimension, 'generators')
    size = 1 << len(generators)
    if size > 64:
        raise ValueError('The phase-mask backend supports at most 64 sign operations')
    context.limit('max_codewords', size)
    group = np.ones((size, dimension), dtype=np.int64)
    for i, row in enumerate(generators):
        group[1 << i:1 << (i+1)] = group[:1 << i]*(1-2*np.array(row, dtype=np.int64))
    return group


def _rows(rows, dimension, weights, norm, context):
    if not isinstance(rows, list) or any(not isinstance(row, list) or len(row) != dimension
            or any(type(x) is not int for x in row) for row in rows):
        raise ValueError('Representatives must be rectangular exact integer rows')
    context.limit('max_codewords', len(rows))
    for row in rows:
        for value in row:
            context.limit('max_rational_bits', abs(value).bit_length())
        if sum(w*x*x for w,x in zip(weights,row)) != norm:
            raise ValueError('Representative has the wrong norm')


def _backend(rows, weights):
    bound = max((abs(x) for row in rows for x in row), default=0)**2*sum(weights)
    # These are integer matrix products. The absolute sum bound covers every
    # multiplication and partial sum; no rounded distance predicate is used.
    if bound < 2**24:
        return np.float32, 'bounded exact binary32'
    if bound < 2**62 and max(weights) < 2**62:
        return np.int64, 'int64 proven overflow-safe'
    return object, 'Python integers'


def phase_kernel(representatives, group, weights, norm, context):
    """Bit m of K[a,b] means <a, g_m b> > norm/2, with signs retained."""
    count = len(representatives)
    context.limit('max_gram_entries', count*count)
    context.limit('max_rank_work', count*count*len(group)*len(weights))
    dtype, backend = _backend(representatives, weights)
    values = np.asarray(representatives, dtype=dtype).reshape(count, len(weights))
    metric = np.asarray(weights, dtype=dtype)
    kernel = np.zeros((count, count), dtype=np.uint64)
    for phase, signs in enumerate(group):
        image = values*np.asarray(signs,dtype=dtype)
        for start in range(0, count, 64):
            dots = (values[start:start+64]*metric) @ image.T
            kernel[start:start+64] |= (dots > norm//2).astype(np.uint64)*np.uint64(1 << phase)
    return kernel, backend


def screen_orbits(rows, weights, norm, representatives, generators, blocker_limit, context):
    """Screen whole catalogue using one source query per distinct orbit.

    The incumbent must be invariant under the supplied sign group. Arbitrary
    sources can use the identity group (an empty generator list). All incumbent
    points are removable, including axes, equators and accessory layers.
    """
    dimension = len(weights)
    integer(blocker_limit, 'blocker_limit', 0)
    _rows(representatives, dimension, weights, norm, context)
    group = sign_group(generators, dimension, context)
    context.limit('max_rank_work', (len(rows)+len(representatives))*len(group)*dimension)
    context.limit('max_gram_entries', len(rows)*len(group))
    known = {tuple(row):i for i,row in enumerate(rows)}
    permutation = np.empty((len(group),len(rows)), dtype=np.int64)
    for phase, signs in enumerate(group):
        for i,row in enumerate(rows):
            image = tuple(x*int(s) for x,s in zip(row,signs))
            if image not in known:
                raise ValueError('Source is not closed under the supplied sign group; use a subgroup or identity')
            permutation[phase,i] = known[image]
    seen, distinct = set(), []
    for row in representatives:
        key = min(map(tuple,(group*np.asarray(row,dtype=object)).tolist()))
        if key not in seen and key not in known:
            seen.add(key)
            distinct.append(list(key))
    context.limit('max_pairs', len(distinct)*len(rows))
    context.limit('max_gram_entries', min(64,len(distinct))*len(rows))
    dtype, backend = _backend(rows+distinct,weights)
    source = np.asarray(rows,dtype=dtype)
    metric = np.asarray(weights,dtype=dtype)
    active, blockers, histogram = [], [], {}
    for start in range(0,len(distinct),64):
        dots = (np.asarray(distinct[start:start+64],dtype=dtype)*metric) @ source.T
        for row,products in zip(distinct[start:start+64],dots):
            blocked = np.flatnonzero(products > norm//2).tolist()
            histogram[len(blocked)] = histogram.get(len(blocked),0)+1
            if len(blocked) <= blocker_limit:
                active.append(row)
                blockers.append(blocked)
    receipt = {'supplied_representatives':len(representatives),'distinct_new_orbits':len(distinct),
               'active_orbits':len(active),'group_size':len(group),
               'blocker_histogram_orbits':{str(k):v for k,v in sorted(histogram.items())},
               'blocker_limit':blocker_limit,'screening_pairs':len(distinct)*len(rows),
               'integer_backend':backend,'catalogue_sha256':digest(representatives)}
    return active, blockers, group, permutation, receipt


def select_partial_orbits(rows, weights, norm, representatives, generators, context, *,
                          blocker_limit=8, candidate_limit=20000, seconds=1, pool_seed=None, on_progress=None):
    """Search finite sign orbits without requiring any orbit to stay intact."""
    integer(candidate_limit,'candidate_limit')
    if pool_seed is not None:
        integer(pool_seed,'pool_seed',0)
    rng = None if pool_seed is None else np.random.default_rng(pool_seed)
    _seconds(seconds,context)
    active, base_blockers, group, permutation, receipt = screen_orbits(
        rows,weights,norm,representatives,generators,blocker_limit,context)
    if on_progress:
        on_progress({'stage':'screened',**receipt})
    # Lowest blocker prices first. The explicit cap is a search restriction,
    # recorded even when it cuts through an orbit. No optimality is inferred.
    order = sorted(range(len(active)),key=lambda i:(len(base_blockers[i]),active[i]))
    candidates, blockers, orbit_ids, phases, used = [], [], [], [], []
    total_eligible, sampled = 0, []
    for orbit in order:
        phase_rows = {}
        for phase,signs in enumerate(group):
            key = tuple(x*int(s) for x,s in zip(active[orbit],signs))
            phase_rows.setdefault(key,phase)
        for point,phase in phase_rows.items():
            total_eligible += 1
            context.limit('max_materialized_points',len(rows)+min(candidate_limit,total_eligible))
            if len(sampled)<candidate_limit:
                sampled.append((orbit,point,phase))
            elif rng is not None:
                j = int(rng.integers(total_eligible))
                if j<candidate_limit:
                    sampled[j]=(orbit,point,phase)
    orbit_lookup = {}
    for orbit,point,phase in sampled:
        if orbit not in orbit_lookup:
            orbit_lookup[orbit]=len(used)
            used.append(orbit)
        candidates.append(list(point))
        blockers.append(permutation[phase,base_blockers[orbit]].tolist())
        orbit_ids.append(orbit_lookup[orbit])
        phases.append(phase)
    count = len(candidates)
    context.limit('max_materialized_points',len(rows)+count)
    context.limit('max_pairs',count*(count-1)//2)
    receipt.update(eligible_points=total_eligible,selected_pool_points=count,
                   candidate_limit=candidate_limit,pool_truncated=count<total_eligible,
                   pool_seed=pool_seed,pool_sampling='lowest blocker price' if rng is None else 'uniform point reservoir after blocker cutoff',
                   all_initial_points_removable=True,whole_orbits_required=False,
                   source_data_sha256=digest({'coordinates':rows,'weights':weights,'norm_squared':norm}))
    selected, removed = [], []
    if candidates:
        kernel, backend = phase_kernel([active[i] for i in used],group,weights,norm,context)
        ids, message = np.array(orbit_ids),np.array(phases,dtype=np.uint64)
        edges, adjacency = [], []
        for start in range(0,count,64):
            relative = message[start:start+64,None] ^ message[None,:]
            bad = (kernel[ids[start:start+64,None],ids[None,:]] >> relative) & np.uint64(1)
            for local,row in enumerate(bad):
                i = start+local
                adjacency.append(int.from_bytes(np.packbits(row.astype(np.uint8),bitorder='little').tobytes(),'little') & ~(1<<i))
                neighbors = np.flatnonzero(row[i+1:])+i+1
                context.limit('max_gram_entries',2*(len(edges)+len(neighbors)))
                edges.extend([i,int(j)] for j in neighbors)
        # Partition each old blocker's candidate neighborhood into actual
        # cliques. sum(x_i in clique) <= deletion_j strengthens the fractional
        # relaxation while preserving every point-level exchange.
        groups = {}
        for i,blocked in enumerate(blockers):
            for j in blocked:
                groups[j] = groups.get(j,0) | (1<<i)
        blocker_cliques = []
        for j,remaining in sorted(groups.items()):
            while remaining:
                pool,clique = remaining,[]
                while pool:
                    bit = pool & -pool
                    i = bit.bit_length()-1
                    clique.append(i)
                    remaining ^= bit
                    pool &= remaining & adjacency[i]
                if len(clique)>1:
                    blocker_cliques.append([j,clique])
        if on_progress:
            on_progress({'stage':'graph_ready',**receipt,'candidate_conflicts':len(edges),
                         'kernel_orbits':len(used),'shared_blocker_cliques':len(blocker_cliques)})
        proposal = blocker_exchange([],{'blockers':blockers,'candidate_conflicts':edges,
                                        'blocker_cliques':blocker_cliques,
                                        'seconds':seconds,'presolve':False},context)
        selected = proposal.data['selected_candidate_indices']
        removed = proposal.data['removed_source_indices']
        receipt.update(candidate_conflicts=len(edges),kernel_orbits=len(used),
                       shared_blocker_cliques=len(blocker_cliques),
                       kernel_integer_backend=backend,kernel_sha256=digest(kernel),
                       solver=proposal.data['solver'],blockers_sha256=digest(blockers),
                       selected_candidate_indices=selected,removed_source_indices=removed)
    removed_set = set(removed)
    receipt.update(gain=len(selected)-len(removed),new_sota=False,
                   scope='Exact output geometry; finite supplied catalogue and explicit cutoffs; numerical search gives no upper bound')
    data = {'coordinates':[row for i,row in enumerate(rows) if i not in removed_set]+[candidates[i] for i in selected],
            'weights':weights,'norm_squared':norm,'search':receipt}
    report = verify_integer_spherical(data,context)
    return Artifact('spherical_code',data,report.details['facts'],report)


def partial_sign_orbits(inputs, params, context):
    check_params(params,('representatives','generators','blocker_limit','candidate_limit','seconds','pool_seed'))
    if len(inputs)!=1 or inputs[0].kind!='spherical_code' or 'coordinates' not in inputs[0].data:
        raise ValueError('partial_sign_orbits requires one explicit integer spherical code')
    source = inputs[0]
    if verify_artifact(source,context).status!='exact':
        raise ValueError('Independently exact source required')
    rows = source.data['coordinates']
    return select_partial_orbits(rows,source.data.get('weights',[1]*len(rows[0])),source.data['norm_squared'],
        params.get('representatives'),params.get('generators',[]),context,
        blocker_limit=params.get('blocker_limit',8),candidate_limit=params.get('candidate_limit',20000),
        seconds=params.get('seconds',1),pool_seed=params.get('pool_seed'))


OPERATORS = {'search.partial_sign_orbits':partial_sign_orbits}
INPUT_KINDS = {'search.partial_sign_orbits':('spherical_code',)}
OUTPUT_KINDS = {'search.partial_sign_orbits':('spherical_code',)}



def rational_image_exchanges(inputs, params, context):
    """Search arbitrary partial unions of exact rational isometric images.

    Each generator is an integer numerator matrix and a positive denominator.
    ``words`` explicitly bounds the group words explored; no finite closure is
    assumed. Every generated point, including every original point, is free
    to be selected or removed by the shared-deletion exchange backend.
    """
    from math import gcd, lcm
    check_params(params, ('generators', 'words', 'blocker_limit', 'candidate_limit', 'seconds', 'pool_seed'))
    if len(inputs) != 1 or inputs[0].kind != 'spherical_code' or 'coordinates' not in inputs[0].data:
        raise ValueError('rational_image_exchanges requires one explicit integer spherical code')
    src = inputs[0]
    if verify_artifact(src, context).status != 'exact':
        raise ValueError('Independently exact source required')
    rows = src.data['coordinates']; dim = len(rows[0]); weights = src.data.get('weights', [1]*dim)
    generators = params.get('generators')
    if not isinstance(generators, list):
        raise ValueError('generators must be a list of rational matrices')
    gen = []
    for item in generators:
        if not isinstance(item, dict) or set(item) != {'numerator', 'denominator'}:
            raise ValueError('each generator requires numerator and denominator')
        m, den = item['numerator'], item['denominator']
        integer(den, 'denominator', 1)
        if not isinstance(m, list) or len(m) != dim or any(not isinstance(r, list) or len(r) != dim or any(type(x) is not int for x in r) for r in m):
            raise ValueError('generator numerator must be a square integer matrix')
        for r in m:
            for x in r: context.limit('max_rational_bits', abs(x).bit_length())
        # A^T W A = den^2 W, evaluated with unbounded Python integers.
        for i in range(dim):
            for j in range(dim):
                if sum(weights[k]*m[k][i]*m[k][j] for k in range(dim)) != (den*den*weights[i] if i == j else 0):
                    raise ValueError('generator is not an isometry for the source metric')
        gen.append((m,den))
    words = params.get('words', [[i] for i in range(len(gen))])
    if not isinstance(words, list) or any(not isinstance(w, list) or any(type(i) is not int or not 0 <= i < len(gen) for i in w) for w in words):
        raise ValueError('words must be lists of generator indices')
    context.limit('max_materialized_points', len(rows)*(len(words)+1))
    context.limit('max_rank_work', dim**3*(len(gen)+sum(map(len,words))) + len(rows)*dim*dim*len(words))
    identity = [[int(i == j) for j in range(dim)] for i in range(dim)]
    transforms = []
    for word in words:
        mat, den = identity, 1
        for index in word:
            other, d = gen[index]
            mat = [[sum(other[i][k]*mat[k][j] for k in range(dim)) for j in range(dim)] for i in range(dim)]
            den *= d
            common = den
            for r in mat:
                for x in r: common = gcd(common, abs(x))
            if common > 1:
                mat = [[x//common for x in r] for r in mat]; den //= common
            context.limit('max_rational_bits', max([den.bit_length()]+[abs(x).bit_length() for r in mat for x in r]))
        transforms.append((mat,den))
    scale = 1
    for _,den in transforms:
        scale = lcm(scale,den)
        context.limit('max_rational_bits', scale.bit_length())
    scaled_source = [[scale*x for x in row] for row in rows]
    for row in scaled_source:
        for value in row: context.limit('max_rational_bits', abs(value).bit_length())
    context.limit('max_rational_bits', (src.data['norm_squared']*scale*scale).bit_length())
    images = []
    for mat,den in transforms:
        factor = scale//den
        for row in rows:
            image = [factor*sum(mat[i][j]*row[j] for j in range(dim)) for i in range(dim)]
            for value in image: context.limit('max_rational_bits', abs(value).bit_length())
            images.append(image)
    result = select_partial_orbits(scaled_source,weights,src.data['norm_squared']*scale*scale,images,[],context,
        blocker_limit=params.get('blocker_limit',len(rows)),candidate_limit=params.get('candidate_limit',20000),
        seconds=params.get('seconds',1),pool_seed=params.get('pool_seed'))
    data = dict(result.data)
    data['search'] = dict(data['search'], generator_count=len(gen), group_words=words,
        common_denominator=scale, construction='individually selectable rational isometric images')
    report = verify_integer_spherical(data,context)
    return Artifact('spherical_code',data,report.details['facts'],report)


OPERATORS['search.rational_image_exchanges'] = rational_image_exchanges
INPUT_KINDS['search.rational_image_exchanges'] = ('spherical_code',)
OUTPUT_KINDS['search.rational_image_exchanges'] = ('spherical_code',)
