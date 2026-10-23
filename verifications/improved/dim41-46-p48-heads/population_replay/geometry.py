#!/usr/bin/env python3
"""Independent exact verification of a 40+3 dimensional kissing construction.

The published input is the pinned extremal even unimodular lattice P48n
(minimum squared norm 6).  Family sizes for first E8 fibres are separately
proved by verify_e8_uniform_masses.py.  This audit reconstructs every
gluing class and every angle inequality, rather than trusting search output.
"""
import argparse
import hashlib
import json
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpz_mat

ROOT = Path(__file__).resolve().parent
OLD = ROOT
SEARCH = ROOT / 'data'
OUT = Path(__file__).resolve().parent


def integral(matrix):
    result = np.array([[int(x) for x in row] for row in matrix.tolist()], dtype=np.int64)
    assert fmpz_mat(result.tolist()) == matrix
    return result


def digest_matrix(matrix):
    return hashlib.sha256(json.dumps(matrix.tolist(), separators=(',', ':')).encode('ascii')).hexdigest()


def geometry():
    p = np.load(OLD / 'data/prep_e8.npz')
    d = np.load(OLD / 'data/e8_complement_reduced.npz')
    h = np.load(SEARCH / 'head_fibres.npz')
    b = d['b'].astype(np.int64)
    q = d['q'].astype(np.int64)
    B = np.array(json.loads((OLD / 'data/ambient_basis.json').read_text())['basis_B'], dtype=np.int64)
    G = p['G'].astype(np.int64)
    assert digest_matrix(B) == '99ae70a671b3b3971c731730e4318f0f80a6052f39ecd7edc821d37087c19b12'
    assert digest_matrix(G) == 'f815a3fc9bd0c2eb46176b4cfdde8160c324028969163628dc02e330cfce2c41'
    assert np.array_equal(B, p['BB'])
    assert np.array_equal(B @ B.T, 96 * G)
    assert fmpz_mat(G.tolist()).det() == 1
    assert np.array_equal(p['T'] @ B, p['B'])
    assert np.array_equal(p['B'] @ p['B'].T, 96 * p['H'])
    assert np.array_equal(4 * p['H'], 3 * p['A'] @ p['A'].T)
    assert abs(fmpz_mat(p['A'].tolist()).det()) == 256
    assert np.all(p['A'] % 2 == p['A'][:, :1] % 2)
    assert np.all(p['A'].sum(axis=1) % 4 == 0)
    assert abs(fmpz_mat(p['U'].tolist()).det()) == 1
    desired = np.r_[np.eye(8, dtype=np.int64), np.zeros((40, 8), dtype=np.int64)]
    assert np.array_equal(p['U'] @ G @ p['T'].T, desired)
    assert np.array_equal(b @ b.T, 96 * q)
    assert not np.any(b @ p['B'].T)
    assert not np.any((b @ B.T) % 96)
    # Prove this is the FULL integral orthogonal complement, using the
    # unimodular basis U, not merely a sublattice with the right rank.
    coeff = integral(fmpz_mat((b @ B.T // 96).tolist()) * fmpz_mat(G.tolist()).inv()
                     * fmpz_mat(p['U'].tolist()).inv())
    assert not np.any(coeff[:, :8])
    assert abs(fmpz_mat(coeff[:, 8:].tolist()).det()) == 1
    assert fmpz_mat(q.tolist()).det() == 6561
    qi3 = integral(3 * fmpz_mat(q.tolist()).inv())
    M = integral(4 * fmpz_mat(p['A'].tolist()).inv() * fmpz_mat(p['B'].tolist()))
    assert np.array_equal(M @ M.T, 1152 * np.eye(8, dtype=np.int64))
    assert np.array_equal(h['body40_basis'], b)
    assert np.array_equal(h['body40_gram'], q)
    assert np.array_equal(h['E8_integer_orthogonal_basis'], M)

    # Independently enumerate doubled E8 vectors of squared length <=32.
    # For separation we need only a short gluing representative in every
    # class.  No shortest-vector or enumeration-completeness assumption is
    # used in deriving the lower bounds on body distances below.
    vectors = []
    def visit(prefix, remaining, parity):
        if len(prefix) == 8:
            if sum(prefix) % 4 == 0:
                vectors.append(prefix)
            return
        for x in range(-6 + parity, 7, 2):
            if x*x <= remaining:
                visit(prefix + [x], remaining - x*x, parity)
    visit([], 32, 0)
    visit([], 32, 1)
    u = np.array(vectors, dtype=np.int64)
    norms = np.sum(u*u, axis=1)
    assert dict(Counter(map(int, norms))) == {0: 1, 8: 240, 16: 2160, 24: 6720, 32: 17520}
    assert not np.any((u @ p['A'].T) % 4)
    pair = u @ p['A'].T // 4
    lifts = pair @ p['U'][:8] @ B
    dots = lifts @ b.T
    assert not np.any(dots % 96)
    c = dots // 96 @ qi3
    normal3 = 3*lifts - c @ b
    assert not np.any(normal3 @ b.T)
    assert np.array_equal(np.sum(normal3*normal3, axis=1), 72*norms)
    labels, classes = np.unique(c % 3, axis=0, return_inverse=True)
    assert len(labels) == 6561
    assert np.array_equal(labels, h['class_labels'])
    mu12 = np.full(6561, 99999, dtype=np.int64)
    np.minimum.at(mu12, classes, norms)
    witnesses = np.full(6561, -1, dtype=np.int64)
    for i, cl in enumerate(classes):
        if norms[i] == mu12[cl] and witnesses[cl] < 0:
            witnesses[cl] = i
    assert np.all(witnesses >= 0)
    delta12 = 72-mu12
    assert np.array_equal(delta12, h['class_minimum_norm_times12'])
    lookup = {row.astype(np.uint8).tobytes(): i for i, row in enumerate(labels)}
    zero = lookup[bytes(40)]
    assert mu12[zero] == 0
    return dict(p=p, b=b, q=q, B=B, M=M, qi3=qi3, labels=labels,
                lookup=lookup, delta12=delta12, mu12=mu12,
                u=u, uclasses=classes, unorms=norms, witnesses=witnesses, h=h)

