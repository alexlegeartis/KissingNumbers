#!/usr/bin/env python3
"""Seed minimal vectors of P_48 and close them under the catalogue automorphisms.

Stage 1 (seeds): exhaustive sweep of all coefficient vectors of support <= 3 with entries in
{-3..3} and of support 4 with entries in {-1,1}; every one of norm exactly 6 is a minimal
vector.  This costs nothing and already lands in many different automorphism orbits.

Stage 2 (closure): breadth-first closure under the 3 catalogue generators of
(SL2(23)xS3):2 <= Aut(P_48p).  An automorphism acts on integer coefficient ROW vectors by
x -> x A with A G A^T = G, so the norm is preserved exactly.  |group| = 145728, so an orbit
on LINES has at most 72864 elements -- unlike Gamma_72 no single orbit is anywhere near big
enough, which is why many seeds are needed.

Output p48_lines.npy, int8 (N,48), canonicalised lines (first nonzero coordinate positive),
every row of norm exactly 6.

    python gen_lines.py [target]
"""
import sys, os, time, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'p48_gram.npy')):
    raise SystemExit(
        "gen_lines.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of P_48 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/regenerate.py rebuilds and checks the class family from the small files\nthat are.")



def canon(X):
    nz = X != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(X[np.arange(len(X)), first]).astype(X.dtype)
    return X * s[:, None]


def norms(X, G):
    return np.einsum('ij,jk,ik->i', X.astype(np.int64), G, X.astype(np.int64))


def sparse_seeds(G):
    n = 48
    out = []
    # support 1 and 2 and 3, coefficients in -3..3
    for k in (1, 2, 3):
        vals = [v for v in range(-3, 4) if v != 0]
        combos = np.array(list(itertools.product(vals, repeat=k)), dtype=np.int64)
        combos = combos[combos[:, 0] > 0]          # kill the +-x duplication
        supports = np.array(list(itertools.combinations(range(n), k)), dtype=np.int64)
        # batch over supports
        for a in range(0, len(supports), 4000):
            S = supports[a:a + 4000]
            X = np.zeros((len(S) * len(combos), n), dtype=np.int64)
            rows = np.repeat(np.arange(len(S)), len(combos))
            for t in range(k):
                X[np.arange(len(X)), S[rows, t]] = np.tile(combos[:, t], len(S))
            nm = norms(X, G)
            out.append(X[nm == 6])
    # support 4, coefficients +-1
    vals = np.array(list(itertools.product([-1, 1], repeat=4)), dtype=np.int64)
    vals = vals[vals[:, 0] > 0]
    supports = np.array(list(itertools.combinations(range(n), 4)), dtype=np.int64)
    for a in range(0, len(supports), 2000):
        S = supports[a:a + 2000]
        X = np.zeros((len(S) * len(vals), n), dtype=np.int64)
        rows = np.repeat(np.arange(len(S)), len(vals))
        for t in range(4):
            X[np.arange(len(X)), S[rows, t]] = np.tile(vals[:, t], len(S))
        nm = norms(X, G)
        out.append(X[nm == 6])
    X = np.vstack(out)
    return canon(X)


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 3_000_000
    G = np.load(os.path.join(HERE, 'p48_gram.npy')).astype(np.int64)
    GE = np.load(os.path.join(HERE, 'p48_gens.npy')).astype(np.int64)
    GEf = [A.astype(np.float32) for A in GE]
    t0 = time.time()
    S = sparse_seeds(G)
    rng = np.random.default_rng(20260820)
    r = rng.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    h = S @ r
    _, k = np.unique(h, return_index=True)
    S = S[k]
    print("sparse seeds: %d minimal lines of support<=4  (%.0fs)" % (len(S), time.time() - t0))

    X = S.astype(np.int8)
    H = X.astype(np.int64) @ r
    Hs = np.sort(H)
    frontier = X
    while len(X) < target:
        imgs = []
        for A in GEf:
            P = frontier.astype(np.float32) @ A
            imgs.append(np.rint(P).astype(np.int64))
        I = canon(np.vstack(imgs))
        hi = I @ r
        _, k = np.unique(hi, return_index=True)
        I, hi = I[k], hi[k]
        pos = np.clip(np.searchsorted(Hs, hi), 0, len(Hs) - 1)
        fresh = Hs[pos] != hi
        I, hi = I[fresh], hi[fresh]
        if len(I) == 0:
            print("closed at %d lines" % len(X))
            break
        assert (norms(I, G) == 6).all(), "closure broke the norm"
        assert int(np.abs(I).max()) < 127
        frontier = I.astype(np.int8)
        X = np.vstack([X, frontier])
        H = np.concatenate([H, hi])
        Hs = np.sort(H)
        print("  %9d lines (+%d)  %.0fs" % (len(X), len(I), time.time() - t0))
    assert (norms(X, G) == 6).all()
    assert len(np.unique(H)) == len(X)
    np.save(os.path.join(HERE, 'p48_lines.npy'), X)
    print("saved p48_lines.npy %s  max|coord| = %d" % (X.shape, int(np.abs(X).max())))


if __name__ == '__main__':
    main()
