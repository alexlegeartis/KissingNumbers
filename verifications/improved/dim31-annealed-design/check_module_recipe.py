#!/usr/bin/env python3
"""Check the frozen 2240-point base and its exact orthogonal orientation.

This checks a construction ingredient, not the final 31D certificate. Use
verify.py for the final certificate. No floating-point arithmetic is used.
"""
from pathlib import Path
import json
import numpy as np


def main():
    data = Path(__file__).resolve().parent / 'data'
    base = np.load(data / 'base_module2240.npz', allow_pickle=False)['points'].astype(np.int64)
    module = np.load(data / 'module251420.npz', allow_pickle=False)['points'].astype(np.int64)
    with np.load(data / 'orientation251420.npz', allow_pickle=False) as z:
        numerator = z['numerator'].astype(np.int64)
        denominator = int(z['denominator'])
    assert base.shape == module.shape == (2240, 15)
    assert numerator.shape == (15, 15) and denominator == 4
    # Bound operands before int64 matrix arithmetic to rule out overflow.
    assert np.max(np.abs(base)) <= 8 and np.max(np.abs(module)) <= 8
    assert np.max(np.abs(numerator)) <= 4
    assert np.array_equal(numerator.T @ numerator, denominator**2 * np.eye(15, dtype=np.int64))
    assert np.array_equal(base @ numerator, denominator * module)
    assert len(np.unique(base, axis=0)) == 2240
    assert len(np.unique(module, axis=0)) == 2240
    assert np.all(np.sum(base * base, axis=1) == 32)
    assert np.all(np.sum(base[:, 8:]**2, axis=1) == 16)
    gram = base @ base.T
    np.fill_diagonal(gram, -100)
    assert np.max(gram) <= 16
    print(json.dumps({'pass': True, 'points': 2240, 'dimension': 15,
        'norm_squared': 32, 'max_distinct_dot': int(np.max(gram)),
        'orientation_denominator': denominator, 'exact_orientation_match': True}, indent=2))

if __name__ == '__main__':
    main()
