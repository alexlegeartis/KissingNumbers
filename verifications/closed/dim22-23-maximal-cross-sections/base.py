"""Common: Leech minimal vectors (norm^2=32), cross-sections Lambda_23, Lambda_22."""
import numpy as np, os, sys, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from leech import load as leech_load

def leech():
    A = leech_load().astype(np.int32)   # 196560 x 24, norm^2 = 32
    assert A.shape == (196560,24)
    return A

def levels(A, u):
    c = A @ u
    return c

def lam23(A=None):
    """Level-0 cross-section by a minimal vector u.  Returns (vecs, u)."""
    if A is None: A = leech()
    u = A[0].copy()
    c = A @ u
    Z = A[c==0]
    return Z, u

def lam22(A=None):
    """Cross-section by a 60-degree pair (inner product 16 = '2' in norm-4 units)."""
    if A is None: A = leech()
    u = A[0].copy()
    c = A @ u
    cand = np.nonzero(c==16)[0]
    w = A[cand[0]].copy()
    d = A @ w
    Z = A[(c==0)&(d==0)]
    return Z, u, w

if __name__ == "__main__":
    A = leech()
    u = A[0]
    c = A @ u
    print("level distribution:", sorted(collections.Counter(c.tolist()).items()))
    Z,_ = lam23(A)
    print("Lambda_23 cross-section:", Z.shape)
    Z22,u1,u2 = lam22(A)
    print("Lambda_22 cross-section:", Z22.shape, "<u1,u2> =", int(u1@u2))
    # internal inner products
    for name,S in (("L23",Z),("L22",Z22)):
        g = S[:200] @ S.T
        print(name, "ip values seen:", sorted(set(g.ravel().tolist())))
