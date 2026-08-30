#!/usr/bin/env python3
"""Exact verification of the Gamma_72 classes, in integer arithmetic only.

For each saved class C (an (m,72) integer array of coefficient vectors with respect to the
catalogue basis of Gamma_72, Gram matrix G):

  1. every row x satisfies x G x^T = 8  -- it is a minimal vector of Gamma_72;
  2. for every pair x != y in C, |x G y^T| <= 2  -- the class condition (|cos| <= 1/4);
     in particular no pair has |x G y^T| = 8, so no two rows are +-each other;
  3. the classes are pairwise DISJOINT as sets of lines, i.e. no line +-x occurs in two of
     them.

Together with the exact rational pairwise table of derive.py these are the only facts the
cap construction needs.

    python verify_classes.py file1.npy|file.npz [file2 ...]
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

def load(path):
    if path.endswith('.npz'):
        z = np.load(path)
        base = os.path.basename(path)
        # Two layouts.  disjoint_classes.npz stores one (m,72) array a class; the batches
        # materialise.py writes stack a whole batch into `lines` with `offs` delimiting the
        # classes, which is ~5x faster to write and read.  Both are checked the same way.
        if 'lines' in z.files and 'offs' in z.files:
            L, off = z['lines'], z['offs']
            idx = z['index'] if 'index' in z.files else np.arange(len(off) - 1)
            return [('%s:class%d' % (base, int(idx[j])), L[off[j]:off[j + 1]])
                    for j in range(len(off) - 1)]
        return [(base + ':' + k, z[k]) for k in sorted(z.files)]
    return [(os.path.basename(path), np.load(path))]

def canon(V):
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]

def main():
    G = np.load(os.path.join(HERE, '..', 'data', 'gamma72_gram.npy')).astype(np.int64)
    assert (G == G.T).all() and (np.diag(G) == 8).all()
    items = []
    for p in sys.argv[1:]:
        items += load(p if os.path.isabs(p) else
                      (p if os.path.exists(p) else os.path.join(HERE, p)))
    total = 0
    allrows = []
    for name, C in items:
        A = C.astype(np.int64)
        assert A.shape[1] == 72, (name, A.shape)
        nm = np.einsum('ij,jk,ik->i', A, G, A)
        assert (nm == 8).all(), (name, "row of norm != 8", sorted(set(int(v) for v in nm)))
        worst = 0
        step = 3000
        for a in range(0, len(A), step):
            M = A[a:a + step] @ G @ A.T
            for t in range(len(M)):
                M[t, a + t] = 0
            worst = max(worst, int(np.abs(M).max()))
        assert worst <= 2, (name, "class condition violated, max |ip| =", worst)
        print("%-28s %5d lines   all norms 8   max |off-diagonal ip| = %d  OK"
              % (name, len(A), worst))
        total += len(A)
        allrows.append(canon(A))
    if len(allrows) > 1:
        S = np.vstack(allrows)
        u = np.unique(S, axis=0)
        print("\ndisjointness: %d lines in total, %d distinct -> %s"
              % (len(S), len(u), "DISJOINT" if len(u) == len(S) else "*** OVERLAP ***"))
        assert len(u) == len(S)
    print("\nsum of class sizes = %d" % total)

if __name__ == '__main__':
    main()
