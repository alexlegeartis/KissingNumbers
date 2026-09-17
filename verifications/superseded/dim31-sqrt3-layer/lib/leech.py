"""Leech lattice in Cohn's coordinates (norm^2 = 32 for minimal vectors), built from lib/golay.py.

build()      the 196560 minimal vectors, int8, shapes (+-4,+-4,0^22), (+-2^8,0^16) on octads with an
             even number of minus signs, (-+3,+-1^23) with signs from a codeword.
in_lattice() membership test for an integer vector: all coordinates of one parity m; the set of
             coordinates congruent to 2 (m = 0) resp. 1 (m = 1) modulo 4 is a codeword; and the
             coordinate sum is 4m modulo 8.  (Conway-Sloane, SPLAG ch. 4 sec. 11.)
"""
import numpy as np, itertools, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from golay import golay24

_G = None


def code():
    global _G
    if _G is None:
        _G = golay24()
    return _G


def build():
    G = code()
    assert dict(collections.Counter(G.sum(1).tolist())) == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    octads = G[G.sum(1) == 8]
    out = []
    for i, j in itertools.combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                v = np.zeros(24, dtype=np.int8); v[i] = si; v[j] = sj; out.append(v)
    signs = [m for m in range(256) if bin(m).count('1') % 2 == 0]
    for oc in octads:
        pos = np.nonzero(oc)[0]
        for m in signs:
            v = np.zeros(24, dtype=np.int8)
            for t, p in enumerate(pos):
                v[p] = -2 if (m >> t) & 1 else 2
            out.append(v)
    base_all = np.where(G == 1, -1, 1).astype(np.int8)
    for ci in range(4096):
        base = base_all[ci]; c = G[ci]
        for j in range(24):
            v = base.copy()
            v[j] = 3 if c[j] == 1 else -3
            out.append(v)
    A = np.array(out, dtype=np.int8)
    assert A.shape == (196560, 24), A.shape
    n2 = (A.astype(np.int32) ** 2).sum(1)
    assert (n2 == 32).all()
    assert len(set(map(tuple, A.tolist()))) == 196560
    return A


def in_lattice(x):
    """True iff the integer vector x (Cohn scaling) lies in the Leech lattice."""
    x = np.asarray(x, dtype=np.int64)
    if x.shape != (24,):
        return False
    m = int(x[0] % 2)
    if not ((x % 2) == m).all():
        return False
    S = ((x % 4) == (2 if m == 0 else 1)).astype(np.int64)
    if not any((S == c).all() for c in code()):
        return False
    return int(x.sum()) % 8 == (4 * m) % 8
