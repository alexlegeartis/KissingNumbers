#!/usr/bin/env python3
"""Phase 2: turn the compact (automorphism, keep-mask) state into explicit line arrays.

build_classes.py stores each class as an automorphism g of Gamma_72 plus a bit-mask saying
which of the 2106 base lines survived the greedy disjointness filter.  This script expands
that into the actual coefficient vectors

    class i = { canon(x g_i) : x in base class, mask_i[x] }

and writes them in batches of 1000 classes as compressed .npz files
    classes_<first>.npz    lines  int8  (sum_j m_j, 72)   all classes of the batch stacked
                           offs   int64 (nclasses+1,)     class j is lines[offs[j]:offs[j+1]]
                           index  int64 (nclasses,)       global class number of each
    (one deflate stream per file rather than one per class: ~5x faster to write and read)

int8 is exact here: coefficients of Gamma_72 minimal vectors never exceed 18 in absolute
value, which is asserted for every single row written.

    python materialise.py [batch] [first] [last]
"""
import numpy as np, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'data')   # the package ships its inputs in data/
STATE = os.environ.get('KN_STATE', os.path.join(HERE, 'state'))
OUTD = os.environ.get('KN_OUT', os.path.join(HERE, 'classes'))


def canon(V):
    """sign convention: the first non-zero coordinate is positive (same as dim73-80)"""
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]


def batches():
    out = []
    for f in sorted(os.listdir(STATE)):
        if f.startswith('gens_') and f.endswith('.npz'):
            out.append((int(f[5:10]), os.path.join(STATE, f)))
    return sorted(out)


def main():
    os.makedirs(OUTD, exist_ok=True)
    G = np.load(os.path.join(SRC, 'gamma72_gram.npy')).astype(np.int64)
    C = np.load(os.path.join(SRC, 'class_best_2106.npy')).astype(np.int64)
    m = len(C)
    lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hi = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    t0 = time.time()
    grand = 0
    maxabs = 0
    for start, path in batches():
        if start < lo or start >= hi:
            continue
        dst = os.path.join(OUTD, 'classes_%05d.npz' % start)
        z = np.load(path)
        gs, kp = z['g'].astype(np.int64), z['keep']
        if os.path.exists(dst):
            grand += int(z['size'].sum())
            continue
        parts, offs = [], [0]
        for j in range(len(gs)):
            g = gs[j]
            assert (g @ G @ g.T == G).all(), 'stored matrix is not an automorphism'
            X = canon(C @ g)
            mask = np.unpackbits(kp[j])[:m].astype(bool)
            X = X[mask]
            a = int(np.abs(X).max())
            maxabs = max(maxabs, a)
            assert a <= 127
            assert len(X) == int(z['size'][j]), 'size disagrees with the build record'
            parts.append(X.astype(np.int8))
            offs.append(offs[-1] + len(X))
            grand += len(X)
        L = np.vstack(parts)
        np.savez_compressed(dst, lines=L, offs=np.array(offs, dtype=np.int64),
                            index=np.arange(start, start + len(gs), dtype=np.int64))
        print('%s  %d classes  %d lines  running total %d  (%.1fs)'
              % (os.path.basename(dst), len(gs), len(L), grand, time.time() - t0), flush=True)
    print('materialised %d lines, max |coefficient| = %d, %.1fs'
          % (grand, maxabs, time.time() - t0))


if __name__ == '__main__':
    main()
