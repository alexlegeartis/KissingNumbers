#!/usr/bin/env python3
"""Global disjointness of ALL 32000 classes.

verify.log covered 31050 classes (46 198 171 lines).  The family was later resumed to 32000
(47 150 230 lines) and that extension was never re-checked -- yet dimension 95 uses 31072
classes, 497 of which have index >= 31050.  So the extension has to be verified before
dimension 95 can stand.

Every line of every class is hashed with a FRESH random vector, independent of the one the
build used; the key is |x . r| with |r_i| < 2^51 and |x_i| <= 20, so |x.r| < 72*20*2^51 < 2^63
and the arithmetic is exact.  Equal lines therefore always give equal keys, and
distinct == total proves no line occurs in two classes."""
import numpy as np, os, sys, time
HERE=os.path.dirname(os.path.abspath(__file__))
CD=os.path.join(HERE,'classes')   # produced by scripts/materialise.py
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 424242)
r=rng.integers(-(2**51),2**51,size=72,dtype=np.int64)
keys=[]; tot=0; mx=0; t0=time.time()
if not os.path.isdir(CD) or len([f for f in os.listdir(CD) if f.endswith('.npz')]) < 32:
    print(
        "SKIP -- the materialised classes are not present.  "
        "This reads the 32000 materialised classes, which are 1.6 GB and are not "
        "shipped.\nBuild them first -- about an hour, and they are a deterministic "
        "function of three\nsmall files that ARE shipped:\n\n"
        "    python scripts/build_classes.py\n"
        "    python scripts/materialise.py\n\n"
        "For a quick check that the family really is regenerable, and that the sizes it "
        "gives\nare the shipped ones, run scripts/verify_rebuild.py instead: a minute on "
        "a prefix.")
    raise SystemExit(0)
for f in sorted(os.listdir(CD)):
    z=np.load(os.path.join(CD,f))
    L=z['lines'].astype(np.int64)
    mx=max(mx,int(np.abs(L).max()))
    k=np.abs(L@r)
    assert k.min()>0, "degenerate hash"
    keys.append(k); tot+=len(L)
    print("  %s -> %d lines cumulative (%.0fs)"%(f,tot,time.time()-t0),flush=True)
K=np.concatenate(keys)
assert 72*mx*(2**51) < 2**63, "hash could overflow"
u=len(np.unique(K))
print()
print("max |coefficient| over all lines = %d ; hash bound 72*%d*2^51 = %.3e < 2^63"%(mx,mx,72*mx*2.0**51))
print("TOTAL %d lines, %d distinct keys -> %s"%(tot,u,"DISJOINT" if u==tot else "COLLISION"))
sizes=np.sort(np.load(os.path.join(HERE,'..','data','class_sizes_32000.npy')).ravel())[::-1]
print("class_sizes.npy: %d entries summing to %d  (matches: %s)"%(len(sizes),sizes.sum(),sizes.sum()==tot))
assert u==tot and sizes.sum()==tot
print()
print("ALL 32000 CLASSES VERIFIED PAIRWISE DISJOINT")
