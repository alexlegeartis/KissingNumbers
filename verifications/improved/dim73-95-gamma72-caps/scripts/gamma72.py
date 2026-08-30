#!/usr/bin/env python3
"""Utilities for Nebe's extremal even unimodular lattice Gamma_72 in explicit coordinates.

Data source: the Catalogue of Lattices entry
    http://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html   ("Gamma72")
parsed by parse_neb72.py into
    gamma72_gram.npy   72x72 int64 Gram matrix G of a basis, all diagonal entries 8
    gamma72_gens.npy   6x72x72 int64 generators of a subgroup SL2(25) x PSL2(7):2 of Aut

Reference: G. Nebe, "An even unimodular 72-dimensional lattice of minimum 8",
J. reine angew. Math. 673 (2012) 237-247 (arXiv:1008.2862).

Lattice vectors are integer coefficient vectors x in Z^72; <x,y> = x^T G y.
Minimal vectors are those with x^T G x = 8.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'gamma72_gram.npy')):
    raise SystemExit(
        "gamma72.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of Gamma_72 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/verify_rebuild.py rebuilds and checks the 32000-class family from the small files\nthat are.")


def load_gram():
    return np.load(os.path.join(HERE, 'gamma72_gram.npy')).astype(np.int64)

def load_gens():
    return np.load(os.path.join(HERE, 'gamma72_gens.npy')).astype(np.int64)


def norms(X, G):
    """squared norms of the rows of X (integer coefficient vectors)"""
    return np.einsum('ij,jk,ik->i', X, G, X)


def descend(X, G, maxit=400):
    """Greedy single-coordinate norm minimisation, batched over the rows of X.

    For a lattice vector x and a basis vector e_k, the norm change of x -> x + t e_k is
        2 t (Gx)_k + t^2 G_kk       with G_kk = 8 for this Gram matrix.
    The best integer t for a given k is t = -round((Gx)_k / 8) and gives the decrease
        d_k = 2 t (Gx)_k + 8 t^2 (<= 0).
    We repeatedly apply the single best coordinate move until no move decreases the norm.
    Returns the resulting X (a local minimum of the norm under coordinate moves).
    """
    X = X.astype(np.int64).copy()
    Gx = X @ G                                    # (M,72)
    M = X.shape[0]
    rows = np.arange(M)
    for _ in range(maxit):
        t = -np.round(Gx / 8.0).astype(np.int64)  # (M,72) best step per coordinate
        d = 2 * t * Gx + 8 * t * t                # (M,72) norm change, <= 0
        k = np.argmin(d, axis=1)
        best = d[rows, k]
        active = best < 0
        if not active.any():
            break
        tk = t[rows, k]
        idx = np.nonzero(active)[0]
        X[idx, k[idx]] += tk[idx]
        Gx[idx] += tk[idx, None] * G[k[idx], :]
    return X


def sample_minimal(G, n, rng, spread=3, maxit=400, chunk=20000):
    """Return an array of distinct minimal LINE representatives (norm 8), roughly n of them.

    Random integer starts are pushed downhill by `descend`; whatever lands on norm 8 is
    kept.  Each line +-x is canonicalised by the sign of its first non-zero coordinate.
    """
    out = {}
    while len(out) < n:
        X = rng.integers(-spread, spread + 1, size=(chunk, 72))
        X = descend(X, G, maxit)
        nm = norms(X, G)
        X = X[nm == 8]
        for x in X:
            nz = np.nonzero(x)[0]
            if len(nz) == 0:
                continue
            if x[nz[0]] < 0:
                x = -x
            out[x.tobytes()] = x
        if len(X) == 0:
            break
    A = np.array(list(out.values())[:n], dtype=np.int64)
    return A


if __name__ == '__main__':
    G = load_gram()
    print("Gram 72x72, diagonal", sorted(set(int(v) for v in np.diag(G))),
          " off-diagonal range", int(G.min()), "..", int(G.max()))
    rng = np.random.default_rng(1)
    import time
    t0 = time.time()
    X = rng.integers(-3, 4, size=(20000, 72))
    Y = descend(X, G)
    nm = norms(Y, G)
    print("descent on 20000 random starts in %.1f s" % (time.time() - t0))
    u, c = np.unique(nm, return_counts=True)
    print("local-minimum norms found:", dict(zip([int(a) for a in u], [int(b) for b in c])))
    print("minimum norm seen:", int(nm.min()))
