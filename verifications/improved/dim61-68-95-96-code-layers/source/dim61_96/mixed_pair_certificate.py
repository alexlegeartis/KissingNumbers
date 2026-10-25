"""Independently verify all support pairs and the resulting 95--96D codes.

No search graph, hash, QR48 membership claim, or search score is trusted.
All supplied supports are checked directly in 96 physical coordinates.
"""
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
try:
    from threadpoolctl import threadpool_limits
except ImportError:
    from contextlib import nullcontext
    def threadpool_limits(limits):
        return nullcontext()

from weighted_layers import (all_binary_words, dense96, kernel, matrix, packed,
                             popcounts64, quaternary24, rref, sign24)
from section95_search import best_rs_section

HERE = Path(__file__).resolve().parent


def main():
    start = perf_counter()
    path = HERE / 'mixed_pair_certificate.npz'
    data = dict(np.load(path, allow_pickle=False))
    low, high = data['supports_low'], data['supports_high']
    assert low.dtype == high.dtype == np.uint64 and low.shape == high.shape
    assert np.all(high < (1 << 32))
    assert np.all(popcounts64(low)+popcounts64(high) == 24)
    assert len(np.unique(np.column_stack((low, high)), axis=0)) == len(low)
    incidence = np.column_stack([((low >> np.uint64(c)) & np.uint64(1)) for c in range(64)] +
                                [((high >> np.uint64(c)) & np.uint64(1)) for c in range(32)]).astype(np.float32)
    rows = incidence.reshape(-1, 24, 4).sum(axis=2)
    assert np.all(np.count_nonzero(rows, axis=1) <= 12)
    maximum = 0
    with threadpool_limits(limits=4):
        for offset in range(0, len(low), 512):
            stop = min(offset+512, len(low))
            # Every product is a sum of at most 96 integers 0 or 1. Each
            # intermediate integer is <=96, so binary32 is exact here.
            products = incidence[offset:stop] @ incidence[:stop].T
            products[np.arange(stop-offset), offset+np.arange(stop-offset)] = 0
            maximum = max(maximum, int(products.max()))
            assert np.all(products <= 12), (offset, maximum)
    sign, _ = sign24()
    assert np.array_equal(sign, data['sign_generator'])
    sign_words = all_binary_words(sign)
    assert len(sign_words) == 16384 and min(popcounts64(sign_words[1:])) == 6
    dual = all_binary_words(matrix(kernel(packed(sign), 24), 24))
    dual_weights = popcounts64(dual)
    assert min(dual_weights[1:]) == 8
    octad = int(dual[np.flatnonzero(dual_weights == 8)[0]])
    retained = np.array([0 if j % 2 else np.count_nonzero(
        popcounts64(sign_words & np.uint64(octad if j == 8 else (1 << j)-1)) == j//2)
        for j in range(9)], dtype=np.int64)
    assert retained.tolist() == [16384, 0, 8192, 0, 6144, 0, 5120, 0, 8960]
    # The sign coordinates are freely reordered on each sparse support.
    # For j<8 the dual distance guarantees the listed uniform projections;
    # for j=8 align the intersection with the exhibited dual octad.
    intersection_sizes = popcounts64(low & np.uint64(255))
    histogram = np.bincount(intersection_sizes, minlength=9)
    added95 = int(histogram @ retained)
    outer, _ = quaternary24()  # Exhaustive 4**9-word distance verification.
    dense = dense96(sign, outer)
    assert len(rref(packed(dense), 96)[0]) == 33
    projection = rref([sum(int(row[c]) << c for c in range(8)) for row in dense], 8)[0]
    projected = all_binary_words(matrix(projection, 8))
    assert len(projection) == 7 and np.all(popcounts64(projected) % 2 == 0)
    dense95 = int(np.count_nonzero(popcounts64(projected) == 4)) * (1 << (33-7))
    assert dense95 == 4697620480
    # The transversal has exactly two active section coordinates; equality
    # of their signs retains half of the distance-six sign code.
    rs = best_rs_section()
    terms95 = {'dense_section': dense95, 'transversal24_section': 1 << 31,
               'new_support24_section': added95, 'support6_section': rs['points'],
               'ordinary_poles': 176}
    terms96 = {'dense': 1 << 33, 'transversal24': 1 << 32,
               'new_support24': len(low)*16384, 'support6': 1 << 21, 'poles': 192}
    before = {95: 7080016652, 96: 13268549824}
    dimensions = {str(d): {'bound': sum(terms.values()), 'previous_bound': before[d],
                          'gain': sum(terms.values())-before[d], 'terms': terms}
                  for d, terms in ((95, terms95), (96, terms96))}
    report = {'status': 'Verified: every physical support pair and exact sign/section counts',
              'supports': len(low), 'previous_supports': 23288,
              'additional_supports': len(low)-23288,
              'support_pairs_checked': len(low)*(len(low)-1)//2,
              'maximum_distinct_support_intersection': maximum,
              'maximum_occupied_rows': int(np.count_nonzero(rows, axis=1).max()),
              'section_coordinates': list(range(8)), 'section_intersection_histogram': histogram.tolist(),
              'retained_sign_counts': retained.tolist(), 'dimensions': dimensions,
              'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'seconds': perf_counter()-start,
              'geometry': {
                  'norm_squared': 96,
                  'dense_min_distance': 'min(4*6, 2*12, 24)=24 gives inner product at most 48',
                  'same_sparse_support': '4*(24-2*6)=48',
                  'different_sparse_supports': '4*12=48, every support pair checked directly',
                  'dense_sparse': '2*24=48',
                  'transversal_sparse': 'one cell per occupied row, at most twelve rows',
                  'support6_sparse': '4*2*6=48',
                  'poles_sparse': '2*sqrt(96)<48',
                  'section95': 'Orthogonal identification of the central hyperplane sum(x_0,...,x_7)=0 with R^95'},
              'scope': 'Constructive lower bounds; no unrestricted optimum or literature-priority claim.'}
    (HERE / 'mixed_pair_certificate.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
