"""Independent integer/Fraction replay of the finite exchange exclusions.

NumPy reads the coordinate archives. No optimizer, floating-point geometric
comparison, baseline all-pairs assertion, or search receipt is trusted here.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(pair):
    require(isinstance(pair, list) and len(pair) == 2, 'rational pair')
    require(all(type(x) is int for x in pair) and pair[1] > 0,
            'integer numerator and positive denominator')
    return Fraction(*pair)


def positive_sqrt3(p, q):
    """Decide p + q sqrt(3) > 0, exactly, including equality."""
    if q == 0:
        return p > 0
    if q > 0:
        return p >= 0 or 3*q*q > p*p
    return p > 0 and p*p > 3*q*q


def dot(x, y):
    return sum(int(a)*int(b) for a, b in zip(x, y))


def load_source(path, dimension):
    with np.load(path, allow_pickle=False) as z:
        if dimension == 26:
            require(int(z['rational_denominator']) == 1, 'rational source denominator')
            require(np.all(z['quadratic_denominators'] == 1)
                    and np.all(z['quadratic_radicands'] == 3), 'quadratic source schema')
            A = 2*np.r_[z['rational_points'], z['quadratic_A']].astype(np.int64)
            B = 2*np.r_[np.zeros_like(z['rational_points']), z['quadratic_B']].astype(np.int64)
        else:
            require(int(z['coordinate_denominator']) == 1
                    and int(z['squared_norm']) == 32
                    and int(z['inner_product_limit']) == 16, 'integer source schema')
            A = 2*z['points'].astype(np.int64)
            B = np.zeros_like(A)
    require(A.shape == B.shape and A.ndim == 2 and A.shape[1] == dimension,
            'source shape')
    require(np.all(np.abs(A) <= 16) and np.all(np.abs(B) <= 16),
            'bounded source coefficients')
    require(np.all(np.sum(A*A+3*B*B, axis=1) == 128)
            and np.all(np.sum(A*B, axis=1) == 0), 'exact source norms')
    return A, B


def verify_math(C, A, B, certificate):
    require(C.dtype.kind in 'iu' and C.ndim == 2 and C.shape[1] == A.shape[1],
            'candidate integer coordinates')
    require(np.all(np.abs(C.astype(np.int64)) <= 16), 'bounded candidate coefficients')
    C = C.astype(np.int64)
    n = len(C)
    require(np.all(np.sum(C*C, axis=1) == 128), 'exact candidate norms')
    require(len(np.unique(C, axis=0)) == n, 'candidate distinctness')
    existing = {row.tobytes() for row in A[~np.any(B, axis=1)]}
    require(all(row.tobytes() not in existing for row in C), 'every candidate is new')
    require(certificate['schema_version'] == 1 and certificate['candidates'] == n,
            'certificate identity')

    new_coefficients = [Fraction(0) for _ in range(n)]
    old_coefficients = {}
    total = Fraction(0)
    owner_comparisons = clique_comparisons = 0
    for row in certificate['cliques']:
        owner, members = row['owner'], row['members']
        require(type(owner) is int and -1 <= owner < len(A), 'source owner index')
        require(isinstance(members, list) and len(members) > 0
                and all(type(i) is int and 0 <= i < n for i in members)
                and len(set(members)) == len(members), 'clique member indices')
        weight = rational(row['weight'])
        require(weight >= 0, 'nonnegative clique weight')
        for at, i in enumerate(members):
            if owner >= 0:
                require(positive_sqrt3(dot(C[i], A[owner])-64,
                                       dot(C[i], B[owner])), 'source owner incompatibility')
                owner_comparisons += 1
            for j in members[:at]:
                require(dot(C[i], C[j]) > 64, 'candidate conflict clique')
                clique_comparisons += 1
            new_coefficients[i] += weight
        if owner < 0:
            total += weight  # sum x_i <= 1
        else:
            old_coefficients[owner] = old_coefficients.get(owner, Fraction(0))-weight
            # sum x_i - y_owner <= 0

    for row in certificate['variable_bounds']:
        kind, index = row['kind'], row['index']
        weight = rational(row['weight'])
        require(weight >= 0 and type(index) is int, 'nonnegative variable-bound weight')
        if kind == 'candidate':
            require(0 <= index < n, 'candidate bound index')
            new_coefficients[index] += weight
        else:
            require(kind == 'deletion' and 0 <= index < len(A), 'source deletion bound index')
            old_coefficients[index] = old_coefficients.get(index, Fraction(0))+weight
        total += weight  # each indicator is at most one

    require(all(v >= 1 for v in new_coefficients), 'all candidate indicators covered')
    require(all(v >= -1 for v in old_coefficients.values()), 'source deletion coefficient bound')
    require(total == rational(certificate['exact_gain_upper_bound']), 'exact weighted bound')
    require(0 <= total < 1 and certificate['integer_gain_upper_bound'] == 0,
            'strict bound below one')
    return dict(status='PASS', candidates=n, source_points=len(A),
                exact_gain_upper_bound=[total.numerator, total.denominator],
                integer_gain_upper_bound=0, cliques=len(certificate['cliques']),
                source_owner_products_checked=owner_comparisons,
                candidate_clique_products_checked=clique_comparisons,
                arithmetic='Integer dot products, exact Q(sqrt(3)) signs, and Fraction weights.',
                scope='The explicit finite pool only; all source points may be deleted. '
                      'No assertion of global optimality or complete source geometry.')


def verify_record(record, manifest):
    source = ROOT/'data'/record['source']
    pool = ROOT/'data'/record['pool']
    proof = ROOT/'data'/record['certificate']
    for path in (source, pool, proof):
        relative = path.relative_to(ROOT).as_posix()
        require(hashlib.sha256(path.read_bytes()).hexdigest() == manifest['sha256'][relative],
                'pinned file hash: '+relative)
    A, B = load_source(source, record['dimension'])
    require(len(A) == record['baseline'], 'source cardinality')
    with np.load(pool, allow_pickle=False) as z:
        C = z['candidates']
        if 'reference_pool' in record:
            reference = ROOT/'data'/record['reference_pool']
            relative = reference.relative_to(ROOT).as_posix()
            require(hashlib.sha256(reference.read_bytes()).hexdigest()
                    == manifest['sha256'][relative], 'reference pool hash')
            with np.load(reference, allow_pickle=False) as original:
                require(np.all(original['candidate_B'] == 0)
                        and np.array_equal(C, 2*original['candidate_A']),
                        'identity with the original octad pool')
                require(np.array_equal(z['blocker_indptr'], original['blocker_indptr'])
                        and np.array_equal(z['blocker_indices'], original['blocker_indices']),
                        'original blocker index identity')
    certificate = json.loads(proof.read_text())
    result = verify_math(C, A, B, certificate)
    require(len(C) == record['candidates'], 'record candidate count')
    return dict(name=record['name'], dimension=record['dimension'], **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    manifest = json.loads((ROOT/'manifest.json').read_text())
    records = [verify_record(r, manifest) for r in manifest['records']]
    report = dict(status='PASS', records=records, seconds=time.monotonic()-started,
                  no_new_lower_bound=True)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
