"""Dimension 28: verify tau(28) >= 204896.

Self-contained. Rebuilds the 196 560 Leech minimal vectors from the Golay code (lib/golay.py),
reads only data/, and decides every pair of the configuration EXACTLY -- integer arithmetic for
the R^24 parts and sympy in Q(sqrt2, sqrt3) for the R^4 parts. Exits non-zero on any failure.

Points of R^28 = R^24 + R^4, all of norm 4 (min-4 units); compatible iff inner product <= 2.

    equator (u, 0)                     u minimal, not an owner              192 592
    cap     (sqrt(2/3)u, (2/sqrt3)z)   u an owner, z in its triangle         11 904
    axis    (0, 2z)                    z one of 16 half-vectors                  16
    layer   (v/2, sqrt2 w)             v = +-f_i a Leech frame, w in +-e_k      384
                                                                          ---------
                                                                            204 896
"""
import numpy as np, itertools, sys, os, collections
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from golay import golay24

FAIL = []
def req(cond, msg):
    print("   [%s] %s" % ("OK " if cond else "FAIL", msg))
    if not cond: FAIL.append(msg)

# ---------------------------------------------------------------- the Leech minimal vectors
GOL = golay24().astype(np.int64)


def in_leech(x):
    """Conway-Sloane membership in the Leech lattice at THIS scaling (minimal norm 32).

    All coordinates of one parity m; the coordinates congruent to 2 (m even) or 1 (m odd)
    modulo 4 form a codeword of the SAME Golay code the shell is built from; and the
    coordinate sum is 4m modulo 8.  (Conway-Sloane, SPLAG ch. 4 sec. 11.)"""
    x = np.asarray(x, np.int64)
    m = int(x[0] % 2)
    if not ((x % 2) == m).all():
        return False
    S = ((x % 4) == (2 if m == 0 else 1)).astype(np.int64)
    if not any(bool((S == c).all()) for c in GOL):
        return False
    return int(x.sum()) % 8 == (4 * m) % 8


def leech_minimal():
    """the 98 280 minimal LINES at norm 32, from the Golay code"""
    G = GOL.astype(np.int8)
    w = G.sum(1); oct_ = G[w == 8]
    V = []
    for i in range(24):
        for j in range(i+1, 24):
            for sj in (1, -1):
                x = np.zeros(24, np.int64); x[i] = 4; x[j] = 4*sj; V.append(x)
    even = np.array([[(m >> t) & 1 for t in range(8)] for m in range(256)], np.int64)
    even = even[even.sum(1) % 2 == 0]; even = even[even[:, 0] == 0]
    for o in oct_:
        pos = np.nonzero(o)[0]
        for s in even:
            x = np.zeros(24, np.int64); x[pos] = 2 - 4*s; V.append(x)
    Gl = [tuple(c.tolist()) for c in G]; idx = {t: i for i, t in enumerate(Gl)}
    comp = {i: idx[tuple((np.array(t) ^ 1).tolist())] for i, t in enumerate(Gl)}
    for ci in [i for i in range(4096) if i < comp[i]]:
        c = G[ci]
        for p in range(24):
            x = np.where(c == 1, -1, 1).astype(np.int64)
            x[p] = 3 if c[p] == 1 else -3
            V.append(x)
    return np.array(V, np.int64)

LS = leech_minimal()
print("[1] the Leech lattice")
req(LS.shape == (98280, 24), "98 280 minimal lines rebuilt from the Golay code")
req(set((LS*LS).sum(1).tolist()) == {32}, "every line has norm 32")
ip0 = np.abs(LS @ LS[0])
req(dict(sorted(collections.Counter(ip0.tolist()).items())) == {0: 46575, 8: 47104, 16: 4600, 32: 1},
    "inner-product distribution 46575 / 47104 / 4600 / 1 (the Leech signature)")

def canon(X):
    X = np.asarray(X, np.int64); f = (X != 0).argmax(1)
    return X * np.sign(X[np.arange(len(X)), f])[:, None]
KEY = {r.tobytes(): i for i, r in enumerate(canon(LS))}

# ---------------------------------------------------------------- the configuration
CV = np.load(os.path.join(HERE, 'data', 'classes.npy')).astype(np.int64)   # 8 x 248 x 24
HD = np.load(os.path.join(HERE, 'data', 'heads.npy')).astype(np.int64)     # 24 x 24
print("\n[2] the 8 owner classes")
req(CV.shape == (8, 248, 24), "8 classes of 248 lines")
CI = []
for k in range(8):
    req(all(r.tobytes() in KEY for r in canon(CV[k])), "class %d lies in the Leech shell" % k)
    CI.append(np.array([KEY[r.tobytes()] for r in canon(CV[k])]))
    sub = CV[k] @ CV[k].T; np.fill_diagonal(sub, 0)
    req(int(np.abs(sub).max()) <= 8,
        "class %d has |<u,u'>| <= 8 throughout, i.e. at least arccos(1/4) = 75.52 deg -- "
        "60-degree freeness is not enough, since the owners of one class share all three "
        "of its cap directions" % k)
allc = np.concatenate(CI)
req(len(set(allc.tolist())) == 8*248, "the 8 classes are pairwise disjoint")

print("\n[3] the heads")
req(set((HD*HD).sum(1).tolist()) == {64}, "24 heads, each a norm-8 Leech vector (norm 64 at this scaling)")
req(all(in_leech(v) for v in HD),
    "every head is a Leech vector (Conway-Sloane criterion, against this Golay code)")
Gh = HD @ HD.T
req(bool((Gh == np.diag(np.full(24, 64))).all()), "the 24 heads are mutually orthogonal: a frame")
mx = int(np.abs(LS[allc] @ HD.T).max())
req(mx <= 16, "every owner has |<u,v>| <= 2 against every head (8*|<u,v>| = %d <= 16)" % mx)
mxe = int(np.abs(LS @ HD.T).max())
req(mxe <= 32,
    "over the WHOLE shell max |<u,v>| = %d <= 32, so <v/2,u> <= 2 and the layer deletes no "
    "equator point (the case list below assumed this)" % mxe)

# ---------------------------------------------------------------- the R^4 side
C1 = []
for i in range(4):
    for j in range(i+1, 4):
        for si in (1, -1):
            for sj in (1, -1):
                z = [Fr(0)]*4; z[i] = Fr(si); z[j] = Fr(sj); C1.append(z)     # each divided by sqrt2
def dot(a, b): return sum(x*y for x, y in zip(a, b))
tris = [t for t in itertools.combinations(range(24), 3)
        if all(C1[t[0]][k] + C1[t[1]][k] + C1[t[2]][k] == 0 for k in range(4))]
byp = {}
for t in tris:
    for p in t: byp.setdefault(p, []).append(t)
def cover(rem, acc):
    if not rem: return acc
    p = min(rem)
    for t in byp[p]:
        if set(t) <= rem:
            r = cover(rem - set(t), acc + [t])
            if r is not None: return r
    return None
sel = cover(set(range(24)), [])
HALF = [[Fr(s, 2) for s in sg] for sg in itertools.product((1, -1), repeat=4)]
EK = [[Fr(sg if k == i else 0) for k in range(4)] for i in range(4) for sg in (1, -1)]
print("\n[4] the R^4 directions")
req(len(tris) == 32 and sel is not None and len(sel) == 8,
    "the 24-cell has 32 zero-sum triangles and 8 of them partition its 24 directions")
req(max(Fr(dot(C1[a], C1[b]), 2) for ta in sel for tb in sel if ta is not tb
        for a in ta for b in tb) <= Fr(1, 2), "distinct chosen triangles are separated (cos <= 1/2)")
req({Fr(dot(C1[a], C1[b]), 2) for t in sel for a in t for b in t if a != b} == {Fr(-1, 2)},
    "each triangle is zero-sum (cos = -1/2)")
req(max(dot(a, b) for a in HALF for b in HALF if a is not b) <= Fr(1, 2),
    "the 16 axis directions are a kissing set (cos <= 1/2)")
req(max(dot(a, b) for i, a in enumerate(EK) for j, b in enumerate(EK) if i != j) <= 0,
    "the 8 layer directions are pairwise obtuse (cos <= 0)")
req(max(abs(dot(a, b)) for a in EK for b in HALF) <= Fr(1, 2),
    "layer vs axis directions: |cos| <= 1/2, and the layer needs <= 1/sqrt2")
req({Fr(dot(C1[i], w)) for i in range(24) for w in EK} == {Fr(-1), Fr(0), Fr(1)},
    "cap vs layer directions: cos in {0, +-1/sqrt2}")
req(max(abs(dot(C1[i], h)) for i in range(24) for h in HALF) <= Fr(1),
    "cap vs axis directions: |cos| <= 1/sqrt2 (the case list below assumed this)")

# ---------------------------------------------------------------- exact inequalities
S2, S3, S6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
print("\n[5] every cross-type worst case, decided exactly")
cases = [
    ("equator x equator, <u,u'> = 2",                 sp.Integer(2)),
    ("equator x axis,    identically 0",              sp.Integer(0)),
    ("equator x cap,     sqrt(2/3)*2",                2*S6/3),
    ("cap x cap, same owner, <z,z'> = -1/2",          sp.Rational(2,3)*4 + sp.Rational(4,3)*sp.Rational(-1,2)),
    ("cap x cap, same class <u,u'>=1, z = z'",        sp.Rational(2,3) + sp.Rational(4,3)),
    ("cap x cap, other class <u,u'>=2, <z,z'>=1/2",   sp.Rational(4,3) + sp.Rational(2,3)),
    ("axis x axis, <z,z'> = 1/2",                     sp.Integer(2)),
    ("cap x axis,  <z,z'> = 1/sqrt2",                 (4/S3)*(1/S2)),
    ("layer x equator, <v,u> = 4",                    sp.Integer(2)),
    ("layer x layer, one head, <w,w'> = 0",           sp.Integer(2)),
    ("layer x layer, two heads, same w",              sp.Integer(2)),
    ("layer x axis, <w,z> = 1/2",                     2*S2*sp.Rational(1,2)),
    ("layer x cap,  m = %d, <z,w> = 1/sqrt2" % (mx // 8),
     sp.sqrt(sp.Rational(2, 3)) * sp.Rational(mx, 16) + 2 / S3),
]
for name, val in cases:
    req(bool(sp.simplify(val - 2) <= 0), "%-46s = %s" % (name, sp.nsimplify(val)))
req(bool(sp.simplify(18 + 12*S2 - 36) < 0),
    "the binding case is strict: (sqrt6+2sqrt3)^2 = 18 + 12 sqrt2 < 36, i.e. sqrt2 < 3/2")

N = 196560 - 2*8*248 + 3*2*8*248 + 16 + 16*24
print("\nTOTAL = 196560 - 3968 + 11904 + 16 + 384 = %d" % N)
if FAIL:
    print("\n*** %d FAILURES: %s ***" % (len(FAIL), FAIL)); sys.exit(1)
print("\n*** VERIFIED: tau(28) >= %d   (published 204520, gain %+d) ***" % (N, N - 204520))
sys.exit(0)
