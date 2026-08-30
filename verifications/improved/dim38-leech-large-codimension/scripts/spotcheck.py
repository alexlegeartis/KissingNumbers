#!/usr/bin/env python3
"""
Independent numerical spot-check of the dimension-38 configuration.

`verify.py` proves the theorem by a case analysis that reduces every pair to an
integer inequality.  This script does something deliberately different and dumber:
it materialises actual 38-dimensional unit vectors for a large random sample of
the configuration -- deliberately over-sampling the pair types that the case
analysis treats separately -- and checks every pairwise inner product numerically.

It is a guard against a conceptual error in the case analysis, not a substitute
for it.

    python scripts/spotcheck.py [n_sample]
"""
import os, sys, itertools, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

def golay_from_basis():
    B = np.array([[int(ch) for ch in line.strip()]
                  for line in open(os.path.join(DATA, "golay_basis.txt")) if line.strip()],
                 dtype=np.uint8)
    words = np.zeros((4096, 24), dtype=np.uint8)
    for m in range(4096):
        bits = np.array([(m >> t) & 1 for t in range(12)], dtype=np.uint8)
        words[m] = (bits @ B) % 2
    return words

def leech_min_vectors(G):
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
            v = base.copy(); v[j] = 3 if c[j] == 1 else -3
            out.append(v)
    return np.array(out, dtype=np.int8)

def canonical_lines(A):
    W = A.astype(np.int32).copy()
    first = np.argmax(W != 0, axis=1)
    sgn = W[np.arange(len(W)), first]
    W[sgn < 0] *= -1
    return np.unique(W, axis=0)

NS = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
rng = np.random.default_rng(0)

G = golay_from_basis()
A = leech_min_vectors(G).astype(np.int64)
L = canonical_lines(A).astype(np.int64)
P = np.loadtxt(os.path.join(DATA, "dim14_directions.txt"), delimiter=",", dtype=np.int64)
T = np.loadtxt(os.path.join(DATA, "triples.txt"), dtype=np.int64)
lab = np.loadtxt(os.path.join(DATA, "class_labels.txt"), dtype=np.int64)
Nrot = np.loadtxt(os.path.join(DATA, "axis_rotation.txt"), dtype=np.int64)
Drot = int(open(os.path.join(DATA, "axis_rotation_D.txt")).read().strip())
keep = np.loadtxt(os.path.join(DATA, "axis_keep.txt"), dtype=np.int64)

K = int(lab.max()) + 1
sizes = np.bincount(lab, minlength=K)
order = np.argsort(-sizes)
use = order[:len(T)]                      # class c -> triple T[i] where use[i] = c

S23 = np.sqrt(2.0 / 3.0)
S13 = np.sqrt(1.0 / 3.0)
NL = np.sqrt(32.0)
ND = np.sqrt(8.0)

pts = []
kinds = []

def cap(u, w):
    v = np.zeros(38)
    v[:24] = S23 * u / NL
    v[24:] = S13 * w / ND
    return v

def axis(a):
    v = np.zeros(38)
    v[24:] = a / np.linalg.norm(a)
    return v

# --- deliberately include every pair type -------------------------------------
# (a) many caps from a handful of triples, so same-direction and same-triple
#     pairs are heavily represented
for i in rng.choice(len(use), size=8, replace=False):
    c = int(use[i]); a, b, d = T[i]
    idx = np.nonzero(lab == c)[0]
    if len(idx) > 40:
        idx = rng.choice(idx, size=40, replace=False)
    for j in idx:
        for w in (P[a], P[b], P[d]):
            for s in (1, -1):
                pts.append(cap(s * L[j], w)); kinds.append("cap")
# (b) a uniform random sample of caps from all over the configuration
for _ in range(NS):
    i = int(rng.integers(len(use))); c = int(use[i])
    idx = np.nonzero(lab == c)[0]
    j = int(idx[rng.integers(len(idx))])
    w = P[T[i][rng.integers(3)]]
    s = 1 if rng.random() < .5 else -1
    pts.append(cap(s * L[j], w)); kinds.append("cap")
# (c) all axis points
Arot = (P @ Nrot.T)[keep] / Drot
for a in Arot:
    pts.append(axis(a)); kinds.append("axis")

X = np.array(pts)
_, uix = np.unique(np.round(X, 9), axis=0, return_index=True)
uix = np.sort(uix)
kinds = [kinds[i] for i in uix]
X = X[uix]                       # the sampler may draw the same point twice
print("sampled %d distinct points (%s)" % (len(X), dict(collections.Counter(kinds))))
nrm = np.linalg.norm(X, axis=1)
print("norms in [%.15f, %.15f]" % (nrm.min(), nrm.max()))
assert abs(nrm.min() - 1) < 1e-12 and abs(nrm.max() - 1) < 1e-12, "not unit vectors"

TOL = 1e-10
worst = -2.0
bad = 0
CH = 2000
for s in range(0, len(X), CH):
    Gm = X[s:s + CH] @ X.T
    for r in range(Gm.shape[0]):
        Gm[r, s + r] = -2.0                       # ignore the diagonal
    w = Gm.max()
    if w > worst:
        worst = w
    bad += int((Gm > 0.5 + TOL).sum())
print("largest inner product between two distinct sampled points: %.15f" % worst)
print("pairs exceeding 1/2 (tolerance %g): %d" % (TOL, bad))
# distinctness
ok = (bad == 0) and (worst <= 0.5 + TOL)
print("=" * 60)
print("SPOT CHECK PASSED" if ok else "SPOT CHECK FAILED")
print("=" * 60)
sys.exit(0 if ok else 1)
