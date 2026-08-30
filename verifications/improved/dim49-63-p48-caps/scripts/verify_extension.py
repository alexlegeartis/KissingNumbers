#!/usr/bin/env python3
"""SUPERSEDED -- kept for the record.  Independent verification of the 130 new P_48
classes over a 1170-class base, i.e. of a 1300-class family.

The family shipped now is 1400, because dimension 63 reads lambda(15) = 1282 classes;
regenerate.py builds that many and verifies every property the claim rests on,
including that the size vector it derives equals the shipped data/class_sizes_p48.npy.
Nothing writes the data/extra_classes_p48.npz this script reads any more, so it
cannot run -- it exits with the message below.  The numbers in it describe the
intermediate stage, not the current family.

What it checked, at that stage:

Checked here without trusting the builder:
  (a) every row of every new class has x G x^T = 6;
  (b) every class is internally valid: |x G y^T| <= 2 for every pair, over the FULL Gram;
  (c) the new classes are pairwise disjoint FROM EACH OTHER;
  (d) the new classes are disjoint FROM THE EXISTING 1170 -- checked directly on canonical
      coordinate tuples, not through the build's hash;
  (e) the extended family's size table is consistent."""
import numpy as np, os, time
import os
_D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
if not os.path.exists(os.path.join(_D, 'disjoint_classes_p48.npz')):
    raise SystemExit(
        "SKIP -- SUPERSEDED.  This checked a 1300-class family: a 130-class extension of "
        "a\n1170-class base.  The family shipped now is 1400, because dimension 63 reads "
        "lambda(15)\n= 1282 of them, and nothing writes the data/extra_classes_p48.npz "
        "this reads any more.\n\n"
        "Its work is done by\n\n"
        "    python scripts/regenerate.py\n\n"
        "which rebuilds all 1400 from three small files and two fixed seeds, verifies "
        "every\nproperty the claim rests on -- each class a 60-degree code, all "
        "pairwise disjoint --\nand checks that the size vector it derives is "
        "exactly the shipped\ndata/class_sizes_p48.npy.")
G = np.load(os.path.join(_D, 'p48p_gram.npy')).astype(np.int64)
Z = np.load(os.path.join(_D, 'disjoint_classes_p48.npz'))
E = np.load(os.path.join(_D, 'extra_classes_p48.npz'))
old=sorted(Z.files); newk=sorted(E.files)
print(__doc__)
def canon(r):
    t=tuple(int(x) for x in r)
    return t if next((x for x in t if x),0)>0 else tuple(-x for x in t)
t0=time.time()
def canon_arr(A):
    """sign-normalise each row so that its first non-zero coordinate is positive"""
    A=np.ascontiguousarray(A.astype(np.int8))
    nz=(A!=0)
    first=np.argmax(nz,axis=1)
    sgn=np.sign(A[np.arange(len(A)),first]).astype(np.int8)
    return np.ascontiguousarray(A*sgn[:,None])
OLD=canon_arr(np.concatenate([Z[k] for k in old]))
worst=0; parts=[]; tot=0
for k in newk:
    C=E[k].astype(np.int64)
    CG=C@G
    d=np.einsum('ij,ij->i',CG,C)
    assert (d==6).all(), "%s has a row of wrong norm"%k
    IP=CG@C.T; np.fill_diagonal(IP,0)
    w=int(np.abs(IP).max()); worst=max(worst,w)
    assert w<=2, "%s violates the class condition (|ip| = %d)"%(k,w)
    parts.append(canon_arr(C)); tot+=len(C)
NEW=np.concatenate(parts)
print("(a),(b) %d new classes, %d lines, all norms 6, worst |off-diagonal| = %d (needs <= 2)"%(
    len(newk),tot,worst),flush=True)
ALL=np.concatenate([OLD,NEW])
V=ALL.view([('', ALL.dtype)]*ALL.shape[1]).ravel()
u=len(np.unique(V))
print("(c),(d) %d existing + %d new = %d canonical lines, %d distinct -> %s"%(
    len(OLD),len(NEW),len(ALL),u,"ALL DISJOINT" if u==len(ALL) else "COLLISION"))
assert u==len(ALL)
sz_old=np.array(sorted((len(Z[k]) for k in old),reverse=True),np.int64)
sz_new=np.array(sorted((len(E[k]) for k in newk),reverse=True),np.int64)
ext=np.sort(np.concatenate([sz_old,sz_new]))[::-1]
np.save('dim49-63/class_sizes_p48_familyB_ext.npy',ext)
print("(e) extended family: %d classes, sizes %d..%d, %d lines  (was %d classes, %d lines)"%(
    len(ext),ext.min(),ext.max(),ext.sum(),len(sz_old),sz_old.sum()))
assert ext.sum()==sz_old.sum()+sz_new.sum()
print()
print("ALL CHECKS PASSED  (%.0fs)"%(time.time()-t0))
