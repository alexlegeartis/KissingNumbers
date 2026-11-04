"""Replay the compact 36D extension obstruction using Python integers only.

Run from any working directory: python discovery/benchmarks/replay_extensions.py
This checks a restricted construction family, not a kissing-number upper bound.
"""
from itertools import product
import json
from pathlib import Path


def replay(path=None):
    path = Path(path) if path else Path(__file__).with_name('all_linear_extensions_20261104.json')
    report = json.loads(path.read_text())
    runs = report.get('runs') if isinstance(report, dict) else None
    if not isinstance(runs, list) or len(runs) != 2 or any(not isinstance(r, dict) for r in runs):
        raise ValueError('Exactly two nonempty module proofs are required')
    if any(type(r.get('factor_index')) is not int for r in runs) or {r['factor_index'] for r in runs} != {0, 1}:
        raise ValueError('Proof must cover factor indices 0 and 1 exactly once')
    results = []
    for run in runs:
        generator = run['generator_without_extension']
        k = len(generator)
        if k != 9 or any(not isinstance(row, list) or len(row) != 17 or
                         any(type(x) is not int or not 0 <= x < 4 for x in row) for row in generator):
            raise ValueError('Expected a 9 by 17 matrix over Z4')
        # Full binary rank proves freeness over Z4: a nonzero message cannot
        # be a zero codeword, which would not obstruct minimum distance.
        pivots = {}
        for row in generator:
            mask = sum((x % 2) << j for j, x in enumerate(row))
            while mask:
                pivot = mask.bit_length() - 1
                if pivot not in pivots:
                    pivots[pivot] = mask
                    break
                mask ^= pivots[pivot]
        if len(pivots) != k:
            raise ValueError('Generator must have full binary row rank')
        constraints = run['obstruction_constraints']
        if not isinstance(constraints, list) or not constraints:
            raise ValueError('An exclusion needs nonempty obstruction constraints')
        weight = (0, 1, 2, 1)
        for constraint in constraints:
            coefficients = constraint['message_coefficients']
            if not isinstance(coefficients, list) or len(coefficients) != k or not any(coefficients) or any(type(x) is not int or not 0 <= x < 4 for x in coefficients):
                raise ValueError('Witness must be a nonzero Z4 message')
            word = [sum(coefficients[i] * generator[i][j] for i in range(k)) % 4
                    for j in range(17)]
            base = sum(weight[symbol] for symbol in word)
            if type(constraint['base_homogeneous_weight']) is not int or base != constraint['base_homogeneous_weight']:
                raise ValueError('Incorrect base Lee weight')
            if type(constraint['required_extension_weight']) is not int or constraint['required_extension_weight'] != 9 - base:
                raise ValueError('Incorrect extension distance requirement')
        eliminated = [0] * len(constraints)
        survivors = 0
        for extension in product(range(4), repeat=k):
            for index, constraint in enumerate(constraints):
                symbol = sum(a * b for a, b in zip(extension, constraint['message_coefficients'])) % 4
                if weight[symbol] < constraint['required_extension_weight']:
                    eliminated[index] += 1
                    break
            else:
                survivors += 1
        if survivors != 0 or sum(eliminated) != 4 ** k:
            raise ValueError('Exclusion proof leaves surviving extensions')
        expected = run['independent_integer_replay']['first_failing_constraint_counts']
        if any(type(x) is not int for x in expected) or eliminated != expected:
            raise ValueError('Elimination counts do not match replay')
        results.append({'factor_index': run['factor_index'], 'extensions_checked': 4 ** k,
                        'constraints': len(constraints), 'survivors': survivors})
    return {'status': 'PASS', 'scope': 'two fixed Z4 modules; all linear one-coordinate extensions',
            'runs': results}


if __name__ == '__main__':
    print(json.dumps(replay(), indent=2))
