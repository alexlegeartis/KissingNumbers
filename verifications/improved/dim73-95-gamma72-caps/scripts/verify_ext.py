#!/usr/bin/env python3
"""Full verification of the 950 EXTENSION classes (indices 31050..31999).

verify.log's exact sample covered 250 classes drawn from the first 31050; the extension added
later was never checked class by class, and dimension 95 uses 497 of it.  Here every one of
the 950 is checked completely -- the full m x m Gram, not a sample of pairs:

    every row x satisfies x G x^T = 8   (a minimal vector of Gamma_72)
    every pair x != y satisfies |x G y^T| <= 2   (the class condition, |cos| <= 1/4)"""
import numpy as np, os, time
HERE=os.path.dirname(os.path.abspath(__file__))
G=np.load(os.path.join(HERE,'..','data','gamma72_gram.npy')).astype(np.int64)
CD=os.path.join(HERE,'classes')   # produced by scripts/materialise.py
worst=0; nchk=0; nl=0; t0=time.time()
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
    lines,offs,index=z['lines'],z['offs'],z['index']
    for i,ci in enumerate(index.tolist()):
        if ci<31050: continue
        C=lines[offs[i]:offs[i+1]].astype(np.int64)
        CG=C@G
        d=np.einsum('ij,ij->i',CG,C)
        assert (d==8).all(), "class %d has a row of wrong norm"%ci
        IP=CG@C.T; np.fill_diagonal(IP,0)
        w=int(np.abs(IP).max())
        assert w<=2, "class %d violates the class condition (|ip| = %d)"%(ci,w)
        worst=max(worst,w); nchk+=1; nl+=len(C)
print("%d extension classes (%d lines) fully checked in %.0fs"%(nchk,nl,time.time()-t0))
print("every row of norm 8; worst |off-diagonal inner product| over the whole extension = %d (needs <= 2)"%worst)
assert nchk==950
print()
print("ALL EXTENSION CLASSES VALID")
