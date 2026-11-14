"""Read pinned PackingStar inputs or the current verified integral Q(sqrt(3)) codes.

Searches use physical squared norm 32 and inner-product threshold 16.
The current certificates store twice these coordinates, with denominator 2.
"""
import numpy as np
if not __debug__:raise RuntimeError('Run without -O: exact input checks require assertions.')


def _certificate(z, dimension):
    assert int(z['radicand']) == 3
    assert int(z['squared_norm']) == 128
    assert int(z['inner_product_limit']) == 64
    assert int(z['coordinate_denominator']) == 2
    A, B = z['rational'], z['radical']
    assert np.issubdtype(A.dtype, np.integer) and np.issubdtype(B.dtype, np.integer)
    assert A.ndim == 2 and A.shape == B.shape and A.shape[1] == dimension
    if np.any(A % 2) or np.any(B % 2):
        raise ValueError('This search requires integral physical coordinates.')
    return A.astype(np.int64) // 2, B.astype(np.int64) // 2


def _norms(A, B, dimension):
    assert A.ndim == 2 and A.shape == B.shape and A.shape[1] == dimension
    assert np.all(np.sum(A*A + 3*B*B, axis=1) == 32)
    assert not np.any(np.sum(A*B, axis=1))
    return A, B


def load_source26(path):
    with np.load(path, allow_pickle=False) as z:
        if 'rational' in z:
            A, B = _certificate(z, 26)
        else:
            assert int(z['rational_denominator']) == 1
            assert np.all(z['quadratic_denominators'] == 1)
            assert np.all(z['quadratic_radicands'] == 3)
            A = np.r_[z['rational_points'], z['quadratic_A']].astype(np.int64)
            B = np.r_[np.zeros_like(z['rational_points']), z['quadratic_B']].astype(np.int64)
    return _norms(A, B, 26)


def load_source27(path):
    with np.load(path, allow_pickle=False) as z:
        if 'rational' in z:
            A, B = _certificate(z, 27)
        else:
            assert int(z['coordinate_denominator']) == 1
            assert np.issubdtype(z['points'].dtype, np.integer)
            A = z['points'].astype(np.int64)
            B = np.zeros_like(A)
    A, B = _norms(A, B, 27)
    if np.any(B):
        raise ValueError('The 27D search requires a rational source.')
    return A
