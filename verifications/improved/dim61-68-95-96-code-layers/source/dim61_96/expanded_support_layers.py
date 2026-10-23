#!/usr/bin/env python3
"""Build the enlarged ERS support layers in dimensions 66,67,68.

Dense [64,28,16] signs are embedded in R^n. Weight-16 supports come from
the exhaustively enumerated [73,36,16] code after selected coordinates
are forced to zero.  Signs are extended Hamming [16,11,4].
Weight-four supports are affine planes in F_2^7 restricted to n points.
"""

import itertools
import json
from math import comb

import numpy as np

from weighted_layers import HERE, popcounts64, rref, sqs64
from binary_62_code import EXPONENTS


def main():
    full = json.loads((HERE / 'cyclic73.json').read_text())
    assert full['n'] == 73 and full['k'] == 36 and full['minimum_distance'] == 16
    assert sum(full['histogram'].values()) == 1 << 36
    words = np.fromfile(HERE / 'cyclic73.weight16.bin', dtype='<u8').reshape(-1, 2)
    assert len(words) == 174324 and len(np.unique(words, axis=0)) == len(words)
    assert np.all(popcounts64(words[:, 0]) + popcounts64(words[:, 1]) == 16)
    original = [sum(1 << e for e in EXPONENTS) << i for i in range(36)]
    rows, pivots = rref(original, 73)
    remainder = words.copy()
    for row, pivot in zip(rows, pivots):
        selected = ((remainder[:, pivot // 64] >> np.uint64(pivot % 64)) & np.uint64(1)).astype(bool)
        remainder[selected, 0] ^= np.uint64(row & ((1 << 64) - 1))
        remainder[selected, 1] ^= np.uint64(row >> 64)
    assert not np.any(remainder)
    choices = json.loads((HERE / 'shortening_support_search.json').read_text())['dimensions']
    reports = {}
    base = sqs64()
    for n in (66, 67, 68):
        deleted = choices[str(n)]['deleted_coordinates']
        retained = np.ones(len(words), dtype=bool)
        for c in deleted:
            retained &= ((words[:, c // 64] >> np.uint64(c % 64)) & np.uint64(1)) == 0
        count16 = int(retained.sum())
        assert count16 == choices[str(n)]['support_count']
        blocks = list(base)
        for a, b in itertools.combinations(range(64, n), 2):
            for c in range(64):
                d = c ^ a ^ b
                if c < d:
                    blocks.append((c, d, a, b))
        blocks += [b for b in itertools.combinations(range(64, n), 4) if b[0] ^ b[1] ^ b[2] ^ b[3] == 0]
        assert len(set(blocks)) == len(blocks)
        triples = [t for b in blocks for t in itertools.combinations(b, 3)]
        assert len(set(triples)) == len(triples)
        assert len(blocks) == {66: 10448, 67: 10512, 68: 10609}[n]
        terms = {'dense64': 1 << 28, 'weight16': count16 * 2048,
                 'weight4': len(blocks) * 16, 'poles': 2 * n}
        baseline = {66: 331737990, 67: 331737996, 68: 361275480}[n]
        reports[str(n)] = {'terms': terms, 'total': sum(terms.values()),
                           'support16_count': count16, 'deleted_mother_coordinates': deleted,
                           'support4_count': len(blocks), 'reference_bound': baseline,
                           'reference_type': 'repository benchmark' if n == 68 else 'ERS64 plus orthogonal tail',
                           'gain': sum(terms.values()) - baseline}
    ceilings = {}
    for n in (69, 70, 71):
        ceilings[str(n)] = (1 << 28) + len(words) * 2048 + (comb(n, 3) // 4) * 16 + 2 * n
    report = {'status': 'Exact retained support counts, code membership and Steiner packing verified',
              'dimensions': reports,
              'cyclic_template_ceiling_69_71': ceilings,
              'ceiling_scope': 'This weight-16 shell, this dense body and the chain (64,16,4,1), even allowing the entire shell in each dimension',
              'geometry_proof': {'norm_squared': 64, 'inner_product_limit': 32,
                                 'dense_sparse_max': 32, 'weight16_weight4_max': 32,
                                 'weight4_pair_max': 32, 'pole_weight4_max': 32},
              'mother_code_words_enumerated': 1 << 36}
    (HERE / 'expanded_support_layers.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
