"""Exact Q(sqrt(3)) comparison used by the discovery programs."""
import numpy as np

def above(p,q):
    p=np.asarray(p,np.int64);q=np.asarray(q,np.int64);z=64-p
    return np.where(q>=0,(z<0)|(3*q*q>z*z),(z<0)&(z*z>3*q*q))
