"""Reconstruct the pinned 210108-point 27D certificate exactly."""
import numpy as np
from leech import build
I=[0,2,4,5,6,10,11,23]
J=[i for i in range(24)if i not in I]
L=2*build().astype(np.int16)
L=L[np.lexsort(L.T[::-1])]
z=np.load('d27_certificate.npz')
p=z['added_integer'].astype('i2')
e=z['deleted_leech_indices'];ok=np.ones(len(L),dtype=bool);ok[e]=False
Q=np.vstack((np.pad(L[ok],((0,0),(0,3))),p))
assert len(Q)==210108 and np.all(np.sum(Q.astype('i4')**2,axis=1)==128)
mid=z['midpoint_rows'];ends=z['midpoint_endpoints']
assert np.all(p[mid,:24]*2==L[ends].sum(1))
np.save('Q.npy',Q)
np.save('L.npy',L)
Q.tofile('Q.bin')
print('built',len(Q),len(p),len(e))
