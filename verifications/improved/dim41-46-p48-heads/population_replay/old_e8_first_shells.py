#!/usr/bin/env python3
"""Recover all 972 old D5 first-shell counts from explicit shortest E8-core heads.

This is a ranking calculation. Uniform norm0/2/4/6 E8 fibre masses come from
the prior moment derivation; the varying norm8 masses are recovered exactly
from the stored 30,240 dual-core heads. Any proposed improved construction
will still receive direct point enumeration and an independent audit.
"""
import itertools,json,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from flint import fmpz_mat

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'record/structural43_v3'
def key(u):return tuple(int(x)%6 for x in u)+(int(sum(u))%12,)

def e8():
    rows=[]
    def rec(a,rem,parity):
        if len(a)==8:
            if sum(a)%4==0:rows.append(a)
            return
        for x in range(-5+parity,6,2):
            if x*x<=rem:rec(a+[x],rem-x*x,parity)
    rec([],32,1) # even entries
    rec([],32,0) # odd entries
    return np.array(rows,dtype=np.int64)
