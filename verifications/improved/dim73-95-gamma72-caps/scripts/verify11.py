#!/usr/bin/env python3
"""Exact verification of the k = 11 cap directions (dimension 83).

The 604-point record configuration of R^11 has coordinates (A + B sqrt2)/3 with A, B integer,
so its Gram lands in Z[sqrt 2] and every comparison must be made there.  Writing
<x,y> * 9 = P + Q sqrt2 and the scaled norm N = 36:

  * 60-degree code:  <x,y> <= N/2      i.e.   P - N/2 <= -Q sqrt2 ;
  * zero-sum triple: x + y + z = 0 and, for each pair, Q = 0 and P = -N/2 ;
  * pair part:       <x,y> <= -N/2     i.e.   P + N/2 <= -Q sqrt2 .

CORRECTED 2026-08-30.  The pair-part line above read `-(P + N/2) >= Q sqrt2`, which is the
same inequality with its sign reversed, and the code matched the prose.  169 of the 176
shipped pair parts were at cos = +1/2 or cos = 0 and passed it.  The units 2T + P = 344 were
right -- CP-SAT says 344 is the EXACT optimum over triples and pairs together, so dimension
83's number does not move -- but the certificate was not one, and the check could not fail.

Comparisons L <= R sqrt2 are decided by sign analysis and L^2 vs 2R^2 -- integers only.
This beats Lambda_11's 146 triples (292 units) with 84 triples + 176 pairs (344 units), which
is what gives dimension 83 a construction of its own instead of being carried by 82."""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, '..', 'data', 'rec11_exact.npz'))
A, B, c = Z['A'].astype(np.int64), Z['B'].astype(np.int64), Z['c'].astype(np.int64)
PZ = np.load(os.path.join(HERE, '..', 'data', 'rec11_parts.npz'))
T, PR = PZ['triples'], PZ['pairs']
n = len(A)
print(__doc__)
P = (A * c) @ A.T + 2 * ((B * c) @ B.T)
Q = (A * c) @ B.T + (B * c) @ A.T
N = int(P[0, 0])
assert (np.diag(P) == N).all() and (np.diag(Q) == 0).all(), "norms are not rational and equal"
half = -N // 2
assert 2 * half == -N


def le(L, R):
    """L <= R*sqrt(2), integers"""
    if R >= 0 and L <= 0: return True
    if R < 0 and L > 0:   return False
    return (L * L <= 2 * R * R) if L > 0 else (L * L >= 2 * R * R)


# 60-degree code: <x,y> <= N/2, i.e. P + half <= -Q sqrt2, tested on EVERY pair.  The
# version this replaced looked only at pairs with P + half > 0, which is not a necessary
# condition when Q > 0 -- P = 10, Q = 10 is 24.1 and would have been skipped.  Nothing in
# this configuration triggered it (its largest inner product is exactly N/2), but a filter
# that happens to be safe on today's data is not a check.
_L = P + half
_R = -Q
_ok = ((_R >= 0) & (_L <= 0))
_ok |= (_L > 0) & (_R > 0) & (_L * _L <= 2 * _R * _R)
_ok |= (_L <= 0) & (_R < 0) & (_L * _L >= 2 * _R * _R)
np.fill_diagonal(_ok, True)
bad = int((~_ok).sum())
print("   %d points, scaled norm %d, sqrt2 part nonzero on %d pairs" % (n, N, int((Q != 0).sum())))
print("   60-degree code, exact in Z[sqrt2]: %s" % (bad == 0))
assert bad == 0
seen = set()
for a, b, e in T.tolist():
    assert (A[a] + A[b] + A[e] == 0).all() and (B[a] + B[b] + B[e] == 0).all()
    for u, v in ((a, b), (a, e), (b, e)):
        assert P[u, v] == half and Q[u, v] == 0, "triple pair not at exactly 120 degrees"
    for x in (a, b, e):
        assert x not in seen; seen.add(x)
# A pair part needs cos <= -1/2, i.e. <x,y> <= -N/2, i.e. P - half <= -Q sqrt2.  The version
# this replaced asserted le(-P + half, Q), which is the OPPOSITE inequality, under an `or`
# whose other branch le(P + half, Q) is a condition on the CONJUGATE P - Q sqrt2 and is
# satisfied by everything in this Gram.  Both passed on 169 of the 176 shipped pairs, which
# sat at cos = +1/2 (161 of them) and cos = 0 (8): they were not parts at all, and the check
# said so for three weeks.  A predicate that cannot fail on the data reads as coverage, so
# this one is shown to fail first, on witnesses taken from the same Gram.
_pos = np.argwhere((P == -half) & (Q == 0))                  # cos = +1/2
_zer = np.argwhere((P == 0) & (Q == 0))                      # cos = 0
assert len(_pos) and len(_zer), "no witness in this Gram to test the pair predicate against"
for _w in (_pos[0], _zer[0]):
    assert not le(int(P[_w[0], _w[1]]) - half, -int(Q[_w[0], _w[1]])),         "the pair test accepts cos = %d/%d, so it is not testing anything" % (
            int(P[_w[0], _w[1]]), N)
print("   pair test REJECTS a cos = +1/2 and a cos = 0 witness from this Gram: True")
for a, b in PR.tolist():
    assert le(int(P[a, b]) - half, -int(Q[a, b])), "pair part not at cos <= -1/2"
    for x in (a, b):
        assert x not in seen; seen.add(x)
S = n - 3 * len(T) - 2 * len(PR)
assert S >= 0 and len(seen) == 3 * len(T) + 2 * len(PR)
print("   %d triples + %d pairs, %d singles ; 3T+2P+S = %d = |W| : %s"
      % (len(T), len(PR), S, 3 * len(T) + 2 * len(PR) + S, 3 * len(T) + 2 * len(PR) + S == n))
print("   units 2T+P = %d  (Lambda_11 gives 292)" % (2 * len(T) + len(PR)))

# ---- the axis layer, if one has been built ------------------------------------------------
# Dimension 83 had none: the rational rotations of scripts/axis_rotate.py compare Gram entries
# as integers and this configuration's are in Z[sqrt2].  The ROTATION may still be rational,
# though -- only the configuration is not -- so scripts/axis_sqrt2.py takes Q = N/D with
# N^T C N = D^2 C and the poles land in the same ring as the directions.  With
# 9 D <x_j, N x_i> = p + q sqrt2 and both norms 4, the 30-degree condition
# 4 <z,a>^2 <= 3 |z|^2 |a|^2 is (p + q sqrt2)^2 <= 972 D^2, i.e.
# (972 D^2 - p^2 - 2 q^2) - 2 p q sqrt2 >= 0.
_lay = os.path.join(HERE, '..', 'data', 'poles_sqrt2_11.npz')
if os.path.exists(_lay):
    import json
    _z = np.load(_lay)
    _N, _D, _kp = _z['N'].astype(object), int(_z['D'][0]), _z['keep']
    _C = np.zeros((A.shape[1], A.shape[1]), dtype=object)
    for _t in range(A.shape[1]):
        _C[_t, _t] = int(c[_t])
    assert ((_N.T @ _C @ _N) == (_D * _D) * _C).all(), "N^T C N != D^2 C: not an isometry"
    _Ai = A.astype(object) @ _N.T
    _Bi = B.astype(object) @ _N.T
    _Ac, _Bc = (A * c).astype(object), (B * c).astype(object)
    _P = _Ac @ _Ai[_kp].T + 2 * (_Bc @ _Bi[_kp].T)      # cap j against pole i
    _Q = _Ac @ _Bi[_kp].T + _Bc @ _Ai[_kp].T
    _a = 972 * _D * _D - _P * _P - 2 * _Q * _Q
    _b = -2 * _P * _Q
    # a + b sqrt2 >= 0, by the signs and then a^2 against 2 b^2
    _ok = ((_a >= 0) & (_b >= 0))
    _ok |= (_a >= 0) & (_b < 0) & (_a * _a >= 2 * _b * _b)
    _ok |= (_a < 0) & (_b > 0) & (2 * _b * _b >= _a * _a)
    _bad = int((~_ok).sum())
    print()
    print("   axis layer: %d poles, rational rotation of denominator %d" % (len(_kp), _D))
    print("   30-degree cap condition in Z[sqrt2], every pole against every cap: %s"
          % ("all %d x %d pairs clear" % _ok.shape if _bad == 0 else "%d FAIL" % _bad))
    assert _bad == 0, "%d pole-cap pairs are inside a cap" % _bad
    _rows = {r['k']: r for r in json.load(open(os.path.join(HERE, '..', 'rows81-95.json')))}
    assert _rows[11]['poles'] == len(_kp), \
        "rows81-95.json says %d poles, this layer has %d" % (_rows[11]['poles'], len(_kp))
    print("   %d = tau(11): AT THE CEILING" % len(_kp) if len(_kp) == 604 else
          "   %d of tau(11) = 604" % len(_kp))

print()
print("ALL CHECKS PASSED")
