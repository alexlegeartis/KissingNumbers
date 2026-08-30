import numpy as np, sys
from base import leech
from lpvertex import unitize, maxnorm
A=leech(); u1=A[0].copy()
u2=A[np.nonzero(A@u1==16)[0][0]].copy()
u3=A[np.nonzero((A@u1==16)&(A@u2==16))[0][0]].copy()
S=A[(A@u1==0)&(A@u2==0)&(A@u3==0)]
print("Lambda_21:",S.shape)
Y=unitize(S,21)
b,z=maxnorm(Y,tries=10,seed=11)
print("Lambda_21: max|z| over P = %.9f  -> min-max |cos| = %.9f"%(b,0.5/b))
print("   (layer pool for Lambda_22 needs max|cos| <= sqrt(3/11) = %.9f)"%( (3/11)**0.5 ))
np.save("lpz_Lambda_21.npy", z)
