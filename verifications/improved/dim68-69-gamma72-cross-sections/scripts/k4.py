#!/usr/bin/env python3
"""The FOUR-point moment LP.  Validated on the Leech, where the answers are known:
D_4 gives Lambda_20 = 17400 and A_4 gives 15540 (the A_j chain is the wrong sublattice from
rank 4 on -- see KNOWLEDGE.md section 77)."""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from kpoint_lp import run3   # the one shared engine, ../../../../common/
M1L={0:93150,1:47104,2:4600}
D4=[[4,-2,0,0],[-2,4,-2,-2],[0,-2,4,0],[0,-2,0,4]]
A4=[[4,-2,0,0],[-2,4,-2,0],[0,-2,4,-2],[0,0,-2,4]]
A3=[[4,-2,0],[-2,4,-2],[0,-2,4]]
for G,name,true in ((A3,"Leech k=3 A_3",27720),(D4,"Leech k=4 D_4",17400),(A4,"Leech k=4 A_4",15540)):
    t0=time.time()
    r=run3(24,196560,4,[1,2],M1L,G,tdes=11)
    v=r[2] if r and r[0]=='cert' else None
    print("  %-16s certified %-9s true %-7s %s   (%.0fs, %s cells)"%(
        name,v,true,"EXACT" if v==true else ("valid" if v is not None and v<=true else "FAIL"),
        time.time()-t0,r[3] if r else '-'),flush=True)
