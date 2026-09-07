"""An independent end-to-end check of the whole 25-dimensional configuration.

    python fullcheck.py [samples]

Reads only data/; nothing here shares code with verify25.py.

Builds the 197567 points explicitly in R^25 and checks pairs.  The equator
against itself is the Leech shell and is known; every other block of pairs is
checked exhaustively, and on top of that a random sample of pairs is drawn from
the whole configuration so that nothing is checked only by the argument that
says it need not be checked.
"""
import sys

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lib"))
from leech import build as leech_build          # noqa: E402

S = 8
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000000
Mi = leech_build().astype(np.int64)
Mf = Mi.astype(np.float64) / np.sqrt(8.0)
X = np.load(os.path.join(HERE, "data", "heads_X.npy"))
U = np.load(os.path.join(HERE, "data", "heads_U.npy")).astype(np.int64)
P = np.load(os.path.join(HERE, "data", "extra_P.npy"))
pos = {Mi[i].tobytes(): i for i in range(len(Mi))}
oi = np.array([pos[u.tobytes()] for u in U])
alive = np.ones(len(Mi), bool)
alive[oi] = False

EQ = np.hstack([Mf[alive], np.zeros((int(alive.sum()), 1))])
EX = np.hstack([P, np.zeros((len(P), 1))])
CP = np.hstack([X, np.ones((len(X), 1))])
CM = np.hstack([X, -np.ones((len(X), 1))])
PO = np.array([[0.0] * 24 + [2.0], [0.0] * 24 + [-2.0]])
blocks = [("equator", EQ), ("extra", EX), ("cap+", CP), ("cap-", CM), ("poles", PO)]
A = np.vstack([b for _, b in blocks])
N = len(A)
print("configuration: %s  total %d" % ([(n, len(b)) for n, b in blocks], N), flush=True)
print("all squared lengths = 4: %s  (max error %.2e)"
      % (bool(np.abs((A * A).sum(1) - 4).max() < 1e-12), float(np.abs((A * A).sum(1) - 4).max())),
      flush=True)
print("all points distinct: %s"
      % (len({np.round(r, 9).tobytes() for r in A}) == N), flush=True)

worst = -9.0
bad = 0
names = []
off = 0
for nm, b in blocks:
    names.append((nm, off, off + len(b)))
    off += len(b)


def pairmax(Aa, Bb, same=False):
    m = -9.0
    for s in range(0, len(Aa), 512):
        g = Aa[s:s + 512] @ Bb.T
        if same:
            for k in range(len(g)):
                g[k, s + k] = -9
        m = max(m, float(g.max()))
    return m


print("\nblock by block (each must be <= 2):")
for i, (n1, b1) in enumerate(blocks):
    for j, (n2, b2) in enumerate(blocks):
        if j < i:
            continue
        if n1 == "equator" and n2 == "equator":
            print("  %-8s x %-8s  the Leech minimal shell, known" % (n1, n2))
            continue
        m = pairmax(b1, b2, same=(i == j))
        worst = max(worst, m)
        if m > 2 + 1e-9:
            bad += 1
        print("  %-8s x %-8s  max %.12f  %s" % (n1, n2, m, "OK" if m <= 2 + 1e-9 else "FAIL"))

rng = np.random.default_rng(0)
print("\nrandom sample of %d pairs from the whole configuration:" % NS, flush=True)
mx = -9.0
nb = 0
step = 2000000
for s in range(0, NS, step):
    k = min(step, NS - s)
    a = rng.integers(0, N, size=k)
    b = rng.integers(0, N, size=k)
    ok = a != b
    v = np.einsum('ij,ij->i', A[a[ok]], A[b[ok]])
    mx = max(mx, float(v.max()))
    nb += int((v > 2 + 1e-9).sum())
print("  max over the sample %.12f, violations %d" % (mx, nb))

print("\n%s   tau(25) >= %d"
      % ("FULL CHECK PASSES" if (bad == 0 and nb == 0) else "FULL CHECK FAILS", N))
sys.exit(0 if (bad == 0 and nb == 0) else 1)
