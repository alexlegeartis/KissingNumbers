"""Exact 60-point ceiling for the fixed integer norm16 tail alphabet in R6.

There are 12 coordinate axes and 240 signed tetrads. On the tetrads put an
edge when the product exceeds 8. The graph is 16-regular. Exact integer
annihilation by a polynomial whose least root is -4 proves A >= -4 I;
Hoffman's bound is then 240*4/(16+4)=48. The axes add at most 12.
No floating point arithmetic is used in the proof.
"""
import argparse
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    rows = []
    for support in combinations(range(6),4):
        for signs in product((-2,2),repeat=4):
            x = np.zeros(6,np.int64)
            x[list(support)] = signs
            rows.append(x)
    T = np.array(rows,np.int64)
    assert T.shape == (240,6)
    assert len(np.unique(T,axis=0)) == 240
    assert np.all(np.sum(T*T,axis=1)==16)
    D = T@T.T
    assert set(np.unique(D)).issubset({-16,-12,-8,-4,0,4,8,12,16})
    A = (D>8).astype(np.int64)
    np.fill_diagonal(A,0)
    assert np.array_equal(A,A.T)
    assert np.all(A.sum(axis=1)==16)
    roots = [-4,-2,0,2,4,8,12,16]
    P = np.eye(240,dtype=np.int64)
    intermediate_bounds = []
    for r in roots:
        oldmax = int(np.max(np.abs(P)))
        # Every column of A has 16 ones; bound all integer partial sums.
        bound = (16+abs(r))*oldmax
        assert bound < np.iinfo(np.int64).max
        P = P@A-r*P
        assert int(np.max(np.abs(P))) <= bound
        intermediate_bounds.append(bound)
    assert not np.any(P)
    # A is real symmetric, so polynomial annihilation places every eigenvalue
    # among these eight real roots. In particular, lambda_min(A) >= -4.
    assert 240*4 == 48*(16+4)
    # Exhibit equality: the complements of a matching are three supports;
    # their 48 sign patterns have pairwise products at most 8.
    supports = [{0,1,2,3},{0,1,4,5},{2,3,4,5}]
    mask = np.array([set(np.flatnonzero(t)) in supports for t in T])
    U = T[mask]
    Q = U@U.T
    np.fill_diagonal(Q,-16)
    assert len(U)==48 and Q.max()==8
    axes = np.r_[4*np.eye(6,dtype=np.int64),-4*np.eye(6,dtype=np.int64)]
    assert np.max(axes@T.T)==8
    report = dict(result="PASS", signed_tetrads=240, degree=16,
                  exact_annihilator_roots=roots, annihilator_zero=True,
                  intermediate_integer_bounds=intermediate_bounds,
                  least_eigenvalue_lower_bound=-4, signed_tetrad_upper_bound=48,
                  axes=12, tail_code_upper_bound=60, equality_witness_count=60,
                  full_16_point_fiber_head_upper_bound=960,
                  scope="Fixed252-vector integer tail alphabet; complete16-point core fibers in8D. Partial fibers, new tails, and other radii are not excluded.")
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
