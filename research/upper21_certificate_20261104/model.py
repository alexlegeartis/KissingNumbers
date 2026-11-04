"""Verified 1024-word Cayley graph for the 21D line-profile family."""
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
OLD=HERE.parent/'compute19_21_20261104'/'data'
z=np.load(OLD/'component21.npz')
C=z['component'].astype(np.int32)
N=z['neighbors'].astype(np.int32)
P=N[0].copy()
assert np.array_equal(N,np.arange(1024)[:,None]^P)
assert np.array_equal(C[np.arange(1024)[:,None]^np.arange(1024)],C[:,None]^C)
