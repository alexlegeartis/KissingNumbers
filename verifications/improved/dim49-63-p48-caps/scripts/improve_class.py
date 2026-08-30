#!/usr/bin/env python3
"""2-improvement local search on a P_48 class (independent set in the |ip| = 3 conflict
graph on minimal lines).

The pool is far too large to hold the conflict graph, so we work in a NEIGHBOURHOOD of the
current class: scan the whole pool once, keep every line that conflicts with at most `tmax`
members of the current class, and build the induced conflict subgraph on class + candidates.
All its edges are genuine conflicts and all its vertices genuine minimal lines, so any
independent set found in it is a genuine class.

Moves (standard 2-improvement local search for maximum independent set):
  (0,1)  a candidate with no conflict to S is inserted;
  (1,2)  two non-adjacent candidates whose unique S-conflict is the same vertex u: drop u,
         insert both (net +1);
  plateau  random (1,1) swaps.
After every round the pool is rescanned so the neighbourhood follows the class.

    python improve_class.py [classfile] [key] [tmax] [rounds] [timelimit] [poolfile] [out]
"""
import sys, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'p48_gram.npy')):
    raise SystemExit(
        "improve_class.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of P_48 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/regenerate.py rebuilds and checks the class family from the small files\nthat are.")

THR = 2.5


def scan(X, XG, Cidx, tmax, block=100000, cchunk=512):
    """conflict count of every pool line against the class, with progressive pruning:
    a row whose count already exceeds tmax is dropped before the next chunk of the class,
    which removes ~99.9% of the pool after the first few chunks."""
    Cx = X[Cidx].astype(np.float32)
    N = len(X)
    keep_idx = []
    for a in range(0, N, block):
        P = XG[a:a + block].astype(np.float32)
        alive = np.arange(len(P), dtype=np.int64)
        cnt = np.zeros(len(P), dtype=np.int32)
        for c0 in range(0, len(Cx), cchunk):
            ipc = P[alive] @ Cx[c0:c0 + cchunk].T
            cnt[alive] += (np.abs(ipc) >= THR).sum(axis=1).astype(np.int32)
            del ipc
            alive = alive[cnt[alive] <= tmax]
            if len(alive) == 0:
                break
        keep_idx.append(alive + a)
    return np.concatenate(keep_idx)


def build_graph(X, XG, V, block=600):
    n = len(V)
    A = X[V].astype(np.float32)
    AG = XG[V].astype(np.float32)
    rows, cols = [], []
    for a in range(0, n, block):
        M = np.abs(AG[a:a + block] @ A.T)
        ii, jj = np.nonzero(M >= THR)
        ii = ii + a
        m = ii != jj
        rows.append(ii[m].astype(np.int32)); cols.append(jj[m].astype(np.int32))
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    order = np.argsort(rows, kind='stable')
    rows, cols = rows[order], cols[order]
    ptr = np.searchsorted(rows, np.arange(n + 1))
    return ptr, cols


def local_search(ptr, adj, n, S0, rng, timelimit=300):
    def nb(v):
        return adj[ptr[v]:ptr[v + 1]]
    inS = np.zeros(n, dtype=bool); inS[S0] = True
    t = np.zeros(n, dtype=np.int32)
    for v in S0:
        t[nb(v)] += 1

    def add(v):
        assert not inS[v] and t[v] == 0
        inS[v] = True; t[nb(v)] += 1

    def drop(u):
        assert inS[u]
        inS[u] = False; t[nb(u)] -= 1

    best = int(inS.sum()); bestS = inS.copy()
    t0 = time.time()
    while time.time() - t0 <= timelimit:
        for v in np.nonzero((~inS) & (t == 0))[0]:
            if (not inS[v]) and t[v] == 0:
                add(int(v))
        cur = int(inS.sum())
        if cur > best:
            best, bestS = cur, inS.copy()
        cand = np.nonzero((~inS) & (t == 1))[0]
        improved = False
        if len(cand):
            uniq = np.empty(len(cand), dtype=np.int64)
            for k, v in enumerate(cand):
                z = nb(v); uniq[k] = z[inS[z]][0]
            order = np.argsort(uniq, kind='stable')
            cand, uniq = cand[order], uniq[order]
            bnds = np.searchsorted(uniq, np.unique(uniq))
            bnds = np.append(bnds, len(uniq))
            for b in range(len(bnds) - 1):
                grp = [int(x) for x in cand[bnds[b]:bnds[b + 1]]]
                if len(grp) < 2:
                    continue
                u = int(uniq[bnds[b]])
                grp = [v for v in grp if (not inS[v]) and t[v] == 1 and inS[u]]
                if len(grp) < 2:
                    continue
                found = None
                for i2 in range(len(grp)):
                    v = grp[i2]
                    nv = set(int(x) for x in nb(v))
                    for j2 in range(i2 + 1, len(grp)):
                        w = grp[j2]
                        if w not in nv:
                            found = (v, w); break
                    if found:
                        break
                if found is None:
                    continue
                v, w = found
                zv = nb(v); zw = nb(w)
                if not (inS[u] and t[v] == 1 and t[w] == 1
                        and int(zv[inS[zv]][0]) == u and int(zw[inS[zw]][0]) == u):
                    continue
                drop(u); add(v); add(w)
                improved = True
        if improved:
            cur = int(inS.sum())
            if cur > best:
                best, bestS = cur, inS.copy()
            continue
        cand = np.nonzero((~inS) & (t == 1))[0]
        if len(cand) == 0:
            break
        for _ in range(200):
            v = int(rng.choice(cand))
            z = nb(v); u = int(z[inS[z]][0])
            drop(u); add(v)
            cand = np.nonzero((~inS) & (t == 1))[0]
            if len(cand) == 0:
                break
    S = np.nonzero(bestS)[0]
    for v in S:
        z = nb(int(v))
        assert not bestS[z].any(), "local_search produced a dependent set"
    return S, best


def main():
    cf = sys.argv[1] if len(sys.argv) > 1 else 'classes_seed0.npz'
    key = sys.argv[2] if len(sys.argv) > 2 else 'c0'
    tmax = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    tlim = int(sys.argv[5]) if len(sys.argv) > 5 else 240
    poolfile = sys.argv[6] if len(sys.argv) > 6 else 'p48_pool.npy'
    outf = sys.argv[7] if len(sys.argv) > 7 else 'class_p48.npy'
    G = np.load(os.path.join(HERE, 'p48_gram.npy')).astype(np.int64)
    X = np.load(os.path.join(HERE, poolfile))
    N = len(X)
    XG = np.empty((N, 48), dtype=np.int8)
    Gi = G.astype(np.int32)
    for a in range(0, N, 400000):
        XG[a:a + 400000] = (X[a:a + 400000].astype(np.int32) @ Gi).astype(np.int8)
    p = os.path.join(HERE, cf)
    C = np.load(p) if cf.endswith('.npy') else np.load(p)[key]
    rng = np.random.default_rng(2024)
    r = rng.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    h = X.astype(np.int64) @ r
    o = np.argsort(h); hs = h[o]
    hc = C.astype(np.int64) @ r
    pos = np.searchsorted(hs, hc)
    Cidx = o[pos]
    assert (X[Cidx] == C).all(), "class not found in pool"
    print("pool %d, start class %d lines" % (N, len(Cidx)))
    for rd in range(rounds):
        t0 = time.time()
        K = scan(X, XG, Cidx, tmax)
        V = np.unique(np.concatenate([Cidx, K]))
        ptr, adj = build_graph(X, XG, V)
        posmap = {int(v): i for i, v in enumerate(V)}
        S0 = np.array([posmap[int(c)] for c in Cidx], dtype=np.int64)
        newS, sz = local_search(ptr, adj, len(V), S0, rng, timelimit=tlim)
        Cidx = V[newS]
        A = X[Cidx].astype(np.int64)
        worst = 0
        for a in range(0, len(A), 4000):
            M = A[a:a + 4000] @ G @ A.T
            for t2 in range(len(M)):
                M[t2, a + t2] = 0
            worst = max(worst, int(np.abs(M).max()))
        assert worst <= 2, "not a class!"
        assert (np.einsum('ij,jk,ik->i', A, G, A) == 6).all()
        np.save(os.path.join(HERE, outf), X[Cidx].astype(np.int8))
        print("round %d: subgraph %d vertices, %d edges -> class %d lines (max|ip|=%d) (%.0fs)"
              % (rd, len(V), len(adj) // 2, len(Cidx), worst, time.time() - t0))


if __name__ == '__main__':
    main()
