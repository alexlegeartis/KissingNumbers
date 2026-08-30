#!/usr/bin/env python3
"""Scale the Gamma_72 disjoint-class family from 120 to N (default 31050) classes.

Method (section 51 of the knowledge base).  C is the verified 2106-line class of Gamma_72
(pairwise |<u,u'>| <= 2 under the Gram G, every row of norm 8).  If g is an automorphism of
Gamma_72 -- an integer matrix with g G g^T = G -- then C g is again such a class, because g
preserves norms and inner products exactly.  Products of random words in the six catalogue
generators give random automorphisms; the images are almost disjoint, and we make them
exactly disjoint by greedily deleting from each new image the lines already used.

The whole family is ~65 million lines x 72 coordinates ~ 9 GB, far too much for RAM, so this
phase never materialises the lines.  It stores only

    key(line +-x) = |x . r|      with r a random integer vector, |r_i| <= 2^51

in an open-addressed hash table of 2^27 slots.  The trick that makes this cheap is that the
key of a row of C g can be computed WITHOUT forming C g:

    (C g) . r = C . (g r)

so one 72x72 matvec plus one (2106,72) matvec replaces a (2106,72)x(72,72) matmul.

EXACTNESS.  Coefficients of minimal vectors of Gamma_72 are bounded by 18 in absolute value
(measured over the 3.39M-line pool, and asserted for every image in phase 2), so
|x . r| <= 72 * 18 * 2^51 = 2.9e18 < 2^63 = 9.2e18: the dot product never overflows int64.
key is therefore an odd function of x, i.e. a genuine function of the LINE +-x.  A hash
collision between two different lines needs (x-y).r = 0 or (x+y).r = 0, which for a fixed
non-zero integer vector d has probability <= 1/(2*2^51+1); over all 2.1e15 pairs the expected
number of collisions is below 1.  A collision can only cause a line to be DELETED when it
did not have to be, so the output is disjoint unconditionally; only the sizes could suffer,
by an expected ~1 line in 65 million.

Output (directory state/), written every BATCH classes and resumable:
    gens_<i>.npz    g   int8   (BATCH,72,72)      the automorphism of each class
                    keep uint8 (BATCH,264)        packed bit-mask of surviving base lines
                    size int32 (BATCH,)
Everything is exactly reconstructible from these plus class_best_2106.npy.

    python build_classes.py [N] [seed]
"""
import numpy as np, os, sys, time, json

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'data')   # the package ships its inputs in data/
STATE = os.environ.get('KN_STATE', os.path.join(HERE, 'state'))

TBITS = 27
TSIZE = 1 << TBITS
TMASK = np.int64(TSIZE - 1)
RBITS = 51
BATCH = 1000
WALK = 1                      # steps of the walk (each step = one pre-mixed group element)
MIXPOOL = 2048                # pre-mixed elements, each a random length-30 word
MIXLEN = 30


def lookup(table, keys):
    """Vectorised open-addressing membership test.  keys must be non-zero and unique."""
    n = len(keys)
    found = np.zeros(n, dtype=bool)
    active = np.arange(n)
    pos = keys & TMASK
    while active.size:
        v = table[pos]
        k = keys[active]
        hit = v == k
        if hit.any():
            found[active[hit]] = True
        fin = hit | (v == 0)
        keep = ~fin
        active = active[keep]
        pos = (pos[keep] + 1) & TMASK
    return found


def insert(table, keys):
    """Vectorised open-addressing insert of keys known to be absent and unique."""
    active = np.arange(len(keys))
    pos = keys & TMASK
    while active.size:
        occ = table[pos] != 0
        fidx, fpos = active[~occ], pos[~occ]
        cidx, cpos = active[occ], (pos[occ] + 1) & TMASK
        if fidx.size:
            # several distinct keys may probe to the same free slot; one wins, rest re-probe
            up, first = np.unique(fpos, return_index=True)
            table[up] = keys[fidx[first]]
            lost = np.ones(fidx.size, dtype=bool)
            lost[first] = False
            cidx = np.concatenate([cidx, fidx[lost]])
            cpos = np.concatenate([cpos, (fpos[lost] + 1) & TMASK])
        active, pos = cidx, cpos


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 31050
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20260819
    os.makedirs(STATE, exist_ok=True)

    G = np.load(os.path.join(SRC, 'gamma72_gram.npy')).astype(np.int64)
    GE = np.load(os.path.join(SRC, 'gamma72_gens.npy')).astype(np.int64)
    C = np.load(os.path.join(SRC, 'class_best_2106.npy')).astype(np.int64)
    m = len(C)
    assert (G == G.T).all() and (np.diag(G) == 8).all()
    assert (np.einsum('ij,jk,ik->i', C, G, C) == 8).all()
    A = C @ G @ C.T
    np.fill_diagonal(A, 0)
    assert int(np.abs(A).max()) <= 2
    for g in GE:
        assert (g @ G @ g.T == G).all()
    print('base class %d lines, all norms 8, max |off-diag ip| = %d; 6 generators verified'
          % (m, int(np.abs(A).max())), flush=True)

    rng = np.random.default_rng(seed)
    r = rng.integers(-(2 ** RBITS), 2 ** RBITS, size=72, dtype=np.int64)
    np.save(os.path.join(STATE, 'hash_r.npy'), r)

    table = np.zeros(TSIZE, dtype=np.int64)
    I = np.eye(72, dtype=np.int64)

    # Pre-mix: MIXPOOL well-separated group elements, each a random word of length MIXLEN.
    # Stepping the walk by one of these costs a single 72x72 matmul yet keeps the overlap
    # |C ^ Cg| between consecutive images at its well-mixed value (measured 0.050 lines).
    tmix = time.time()
    MIX = np.empty((MIXPOOL, 72, 72), dtype=np.int64)
    for i in range(MIXPOOL):
        a = I.copy()
        for _ in range(MIXLEN):
            a = a @ GE[rng.integers(0, 6)]
        MIX[i] = a
    print('pre-mixed %d group elements (len %d words) in %.1fs, max|entry| %d'
          % (MIXPOOL, MIXLEN, time.time() - tmix, int(np.abs(MIX).max())), flush=True)

    # ---- resume from whatever batches already exist ---------------------------------
    done = 0
    sigs = set()
    t0 = time.time()
    while os.path.exists(os.path.join(STATE, 'gens_%05d.npz' % done)):
        z = np.load(os.path.join(STATE, 'gens_%05d.npz' % done))
        gs, kp = z['g'].astype(np.int64), z['keep']
        for j in range(len(gs)):
            h = C @ (gs[j] @ r)
            keys = np.abs(h)
            mask = np.unpackbits(kp[j])[:m].astype(bool)
            sigs.add(int(np.sum(keys.astype(np.uint64))))
            sel = np.unique(keys[mask])
            insert(table, sel)
        done += len(gs)
        print('resumed batch -> %d classes (%.1fs)' % (done, time.time() - t0), flush=True)
    if done:
        g = np.load(os.path.join(STATE, 'gens_%05d.npz' % (done - len(gs))))['g'][-1].astype(np.int64)
        rng = np.random.default_rng(seed + done)
    else:
        g = I.copy()

    # ---- main loop -----------------------------------------------------------------
    bg = np.zeros((BATCH, 72, 72), dtype=np.int8)
    bk = np.zeros((BATCH, (m + 7) // 8), dtype=np.uint8)
    bs = np.zeros(BATCH, dtype=np.int32)
    nb = 0
    total = int((table != 0).sum()) if done else 0
    dup_images = 0
    maxg = 0
    t0 = time.time()
    tlast = t0
    while done < N:
        # random walk on Aut(Gamma_72): g <- g * (WALK random generators)
        for _ in range(WALK):
            g = g @ MIX[rng.integers(0, MIXPOOL)]
        h = C @ (g @ r)
        keys = np.abs(h)
        assert keys.min() > 0, 'degenerate hash (some line has x.r == 0)'
        sig = int(np.sum(keys.astype(np.uint64)))
        if sig in sigs:                      # same image set as an earlier class
            dup_images += 1
            continue
        sigs.add(sig)
        uk, ui = np.unique(keys, return_index=True)
        old = lookup(table, uk)
        newk = uk[~old]
        insert(table, newk)
        mask = np.zeros(m, dtype=bool)
        mask[ui[~old]] = True
        gm = int(np.abs(g).max())
        maxg = max(maxg, gm)
        assert gm <= 127, 'automorphism entry too large for int8 storage'
        bg[nb] = g.astype(np.int8)
        bk[nb] = np.packbits(mask)
        bs[nb] = int(mask.sum())
        total += int(mask.sum())
        nb += 1
        done += 1
        if nb == BATCH or done == N:
            np.savez_compressed(os.path.join(STATE, 'gens_%05d.npz' % (done - nb)),
                                g=bg[:nb], keep=bk[:nb], size=bs[:nb])
            now = time.time()
            print('%6d classes  total lines %10d  batch sizes %d..%d mean %.1f  '
                  'dup-images %d  max|g| %d  %.1fs (+%.1fs)'
                  % (done, total, int(bs[:nb].min()), int(bs[:nb].max()),
                     float(bs[:nb].mean()), dup_images, maxg, now - t0, now - tlast),
                  flush=True)
            tlast = now
            nb = 0
    print('DONE: %d classes, %d lines, %d duplicate images rejected, %.1fs'
          % (done, total, dup_images, time.time() - t0), flush=True)
    json.dump({'n': done, 'total_lines': total, 'dup_images': dup_images,
               'seed': seed, 'walk': WALK, 'rbits': RBITS, 'tbits': TBITS,
               'seconds': time.time() - t0},
              open(os.path.join(STATE, 'build_summary.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
