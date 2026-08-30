#!/usr/bin/env python3
"""Is S(lambda) improvable, or is the family already where the method must land?

S(lambda) is the maximum number of lines covered by lambda pairwise disjoint classes.  The
ceiling is lambda * alpha with alpha the largest class (7069 achieved).  At lambda = 1282 the
family reaches 7 457 761 against a ceiling of 9 062 458 -- 82.3%.  Is the missing 17.7% search
slack, or is it structural?

MODEL.  Place classes of size a = 7069 one at a time into the shell of N = 26 208 000 lines.
An unstructured class overlaps the already-used set U in about a*U/N lines, so

    dU/di = a (1 - U/N)   =>   U(i) = N (1 - exp(-a i / N)) .

That is the value ANY family of independently-chosen maximal classes must land on; the ceiling
lambda*a is reached only by a family whose classes are disjoint BY CONSTRUCTION."""
import numpy as np, os
N=26208000; a=7069
S=np.sort(np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                     'class_sizes_p48_familyB_ext.npy')).ravel())[::-1].astype(np.int64)
cum=np.concatenate([[0],np.cumsum(S)])
print(__doc__)
print("  lambda    achieved     collision model    ceiling lambda*a    achieved/model  achieved/ceiling")
for L in (1,120,378,577,966,1282,1300):
    L=min(L,len(S))
    mod=N*(1-np.exp(-a*L/N))
    print("  %6d %11d %18.0f %19d      %6.1f%%          %6.1f%%"%(
        L,int(cum[L]),mod,L*a,100.0*cum[L]/mod,100.0*cum[L]/(L*a)))
print()
print("The family tracks the collision model to within a few percent at every lambda, so the")
print("17.7%% shortfall at lambda = 1282 is NOT search slack -- it is the collision loss that")
print("any family of independently chosen classes must pay.  Closing it needs a family whose")
print("classes are disjoint by construction (an H-free class under a subgroup H of Aut(P_48)),")
print("not a longer search.")
