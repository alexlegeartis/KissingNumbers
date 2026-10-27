"""Exact mixed-tail, weighted-priority existence bounds in 77--79 dimensions.

Each vertex is a (norm-ten head line, oriented tail) and represents an
antipodal pair of points. Independent exponential priorities of rate r(v)
select a vertex first in its closed neighborhood with probability
r(v) / sum(r(w) for w in N[v]). Selected vertices are independent.
Rates below are fixed positive integers; no numerical optimizer is used.
See README.md for the geometric and probabilistic proof.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import sympy as sp

from norm10_bound_certificate import theta_coefficients, solve_moments

HERE = Path(__file__).resolve().parent
SPECS = {
    5: (Q(16, 25), 157451),
    6: (Q(16, 25), 256549),
    7: (Q(12, 19), 373462),
}


def candidate_tails(data, rank):
    """Fixed tails and integer rates, defined without search state."""
    old = data['tails%d' % rank].astype(np.int64)
    poles = data['poles%d' % rank].astype(np.int64)
    if rank == 5:
        # Norm-eight D5 vectors in unscaled integer coordinates; multiplying
        # by two gives the additional norm-32 tails in the 2E8 model.
        extra = np.array([tuple(2*x for x in v) + (0, 0, 0)
                          for v in product(range(-2, 3), repeat=5)
                          if sum(x*x for x in v) == 8 and sum(v) % 2 == 0],
                         dtype=np.int64)
        assert len(extra) == 200
        z = np.vstack((extra, old))
        rates = np.array([1050 if np.count_nonzero(v) == 2 else 818 for v in extra]
                         + [882 if np.max(np.abs(v)) == 6 else 1000 for v in old],
                         dtype=np.int64)
    elif rank == 6:
        z = np.vstack((poles, old))
        rates = np.array([1015]*len(poles) + [1000]*len(old), dtype=np.int64)
    else:
        doubled_roots = {tuple(2*v) for v in poles}
        z = old
        rates = np.array([1000 if tuple(v) in doubled_roots else 311 for v in z],
                         dtype=np.int64)
    return z, rates, poles


def check_geometry(z, rates, poles, rank):
    assert z.ndim == 2 and z.shape[1] == 8
    assert rates.shape == (len(z),) and np.all(rates > 0)
    assert len(np.unique(z, axis=0)) == len(z)
    weighted = {tuple(v): int(r) for v, r in zip(z, rates)}
    assert all(weighted.get(tuple(-v)) == r for v, r in zip(z, rates))
    norms = np.sum(z*z, axis=1)
    assert set(norms.tolist()) == ({32, 40} if rank == 5 else
                                  {8, 40} if rank == 6 else {32})
    # The Gram rank equals the Euclidean span rank, and includes the poles.
    together = np.vstack((z, poles))
    assert sp.Matrix((together.T @ together).tolist()).rank() == rank
    assert len(np.unique(poles, axis=0)) == len(poles)
    assert np.all(np.sum(poles*poles, axis=1) == 8)
    pd = poles @ poles.T
    np.fill_diagonal(pd, 0)
    assert np.max(pd) <= 4
    # These conservative bounds also certify every int64 square and weighted
    # sum used below. The final rational sums use Python integers throughout.
    assert np.max(np.abs(z)) <= 6 and len(z) <= 5000
    assert np.max(rates) <= 1050
    return norms


def closed_profiles(z, rates, norms, scale, shell):
    """Check every tail pair; do not assume transitivity or equal degrees."""
    p, d = scale.numerator, scale.denominator
    residual = 8*d - 10*p
    assert Q(3, 5) < scale <= Q(16, 25)
    assert 3*scale + (8-10*scale) <= 4
    assert 10*scale - (8-10*scale) > 4  # all same-head vertices conflict
    assert all(4*d-j*p > 0 for j in (4, 5, 6))
    total_rate = int(sum(map(int, rates)))
    profiles = Counter()
    for start in range(0, len(z), 128):
        dots = z[start:start+128] @ z.T
        norm_products = norms[start:start+128, None] * norms[None, :]
        closed = np.full(len(dots), total_rate, dtype=np.int64)
        for j in (4, 5, 6):
            # residual*dot/sqrt(qr) > 4*d-j*p. Both sides are positive
            # exactly when dot>0, so squaring here preserves strictness.
            conflicts = ((dots > 0) &
                         (residual**2*dots*dots > (4*d-j*p)**2*norm_products))
            closed += (shell[j]//2) * (conflicts @ rates)
        profiles.update((int(rate), int(degree))
                        for rate, degree in zip(rates[start:start+128], closed))
    return profiles


def certificate():
    theta, coefficients = theta_coefficients()
    n8, n10 = theta[4:6]
    assert (n8, n10) == (6218175600, 15281788354560)
    small = solve_moments(n8, 80, {})
    shell = solve_moments(n10, 100, {6: small[4], 10: 2})
    assert sum(shell.values()) == n10
    assert all(shell[j] % 2 == 0 for j in (4, 5, 6))
    source = HERE / 'norm10_tails.npz'
    rows = {}
    with np.load(source, allow_pickle=False) as data:
        for rank, (scale, expected_pairs) in SPECS.items():
            z, rates, poles = candidate_tails(data, rank)
            norms = check_geometry(z, rates, poles, rank)
            profiles = closed_profiles(z, rates, norms, scale, shell)
            alpha = (n10//2) * sum((Q(count*rate, closed)
                                  for (rate, closed), count in profiles.items()), Q(0))
            pairs = -(-alpha.numerator//alpha.denominator)
            assert pairs == expected_pairs
            rows[str(72+rank)] = {
                'bound': n8 + 2*pairs + len(poles),
                'equator': n8,
                'head_scale_squared': str(scale),
                'tail_vectors': len(z),
                'norm_rate_populations': [list(key)+[count] for key, count in
                                          sorted(Counter(zip(map(int, norms),
                                                             map(int, rates))).items())],
                'guaranteed_antipodal_pairs': pairs,
                'pure_tail_points': len(poles),
                'priority_expectation_numerator': str(alpha.numerator),
                'priority_expectation_denominator': str(alpha.denominator),
                'weighted_closed_profiles': [[rate, closed, count]
                                             for (rate, closed), count in sorted(profiles.items())],
                'tail_pair_products_checked': len(z)**2,
            }
    return {
        'status': 'PASS',
        'certificate_type': 'Exact weighted-priority finite-graph existence proof; no full head list',
        'theta_coefficients': theta,
        'modular_form_coefficients': coefficients,
        'minimal_shell_absolute_products_against_norm10': small,
        'norm10_absolute_product_distribution': shell,
        'dimensions': rows,
        'tails_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'literature_inputs': [
            'Gamma72 even unimodular, minimum 8 (Nebe 2010)',
            'All shells are spherical 11-designs (Venkov; Nebe survey Theorem 4.1)',
            'Classical modular-form description of even unimodular theta series',
        ],
    }


def main():
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    result = certificate()
    path = HERE / 'weighted_norm10_certificate.json'
    if '--check' in sys.argv:
        assert json.loads(path.read_text()) == json.loads(json.dumps(result)), 'stale certificate'
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
