"""Cache-resident pair kernel for independent CPU searches; never a verifier."""
import numpy as np
from numba import njit


@njit(cache=True,fastmath=True)
def candidate_pairs(x,threshold):
    n,d=x.shape
    ii=np.empty(n*(n-1)//2,dtype=np.int32)
    jj=np.empty_like(ii)
    count=0
    for i in range(n):
        for j in range(i):
            dot=0.
            for k in range(d):dot+=x[i,k]*x[j,k]
            if dot>threshold:
                ii[count]=i;jj[count]=j;count+=1
    return ii[:count].copy(),jj[:count].copy()


@njit(cache=True,fastmath=True)
def sparse_energy_gradient(x,ii,jj,target,power):
    n,d=x.shape
    grad=np.zeros_like(x)
    values=np.empty(len(ii))
    peak=0.
    for e in range(len(ii)):
        i,j=ii[e],jj[e]
        dot=0.
        for k in range(d):dot+=x[i,k]*x[j,k]
        v=max(dot-target,0.)
        values[e]=v
        peak=max(peak,v)
    if peak==0.:return 0.,grad
    s=0.
    for e in range(len(ii)):
        v=values[e]
        if v==0.:continue
        i,j=ii[e],jj[e]
        if power==2:
            w=v;s+=v*v
        else:
            z=v/peak;w=z**(power-1);s+=w*z
        for k in range(d):
            grad[i,k]+=w*x[j,k]
            grad[j,k]+=w*x[i,k]
    if power==2:
        grad*=2.
        return s,grad
    grad*=2*peak*s**(2./power-1.)
    return peak*peak*s**(2./power),grad


class NeighborList:
    """Safely omit inactive pairs using |delta dot| <= ||dx_i|| + ||dx_j||.

    Rebuild as soon as twice the maximum displacement reaches the saved gap.
    This bound applies to unit vectors in every ambient dimension used here.
    An extra 1e-10 guard is much larger than floating point dot-product error.
    """
    def __init__(self,target,margin=.025):
        self.target=target
        self.margin=margin
        self.reference=None
        self.rebuilds=0

    def evaluate(self,x,power=2):
        if (self.reference is None or
            2*np.sqrt(np.max(np.sum((x-self.reference)**2,axis=1)))>=self.margin-1e-10):
            self.ii,self.jj=candidate_pairs(x,self.target-self.margin-1e-10)
            self.reference=x.copy()
            self.rebuilds+=1
        return sparse_energy_gradient(x,self.ii,self.jj,self.target,power)


@njit(cache=True, fastmath=True)
def pair_energy_gradient(x, target):
    n,d=x.shape
    grad=np.zeros_like(x)
    energy=0.
    for i in range(n):
        for j in range(i):
            dot=0.
            for k in range(d):
                dot+=x[i,k]*x[j,k]
            v=dot-target
            if v>0:
                energy+=v*v
                twice=2*v
                for k in range(d):
                    grad[i,k]+=twice*x[j,k]
                    grad[j,k]+=twice*x[i,k]
    return energy,grad


@njit(cache=True, fastmath=True)
def power_energy_gradient(x, target, power):
    n,d=x.shape
    grad=np.zeros_like(x)
    # Relative scaling keeps high powers numerically well conditioned.
    peak=0.
    for i in range(n):
        for j in range(i):
            dot=0.
            for k in range(d):
                dot+=x[i,k]*x[j,k]
            peak=max(peak,dot-target)
    if peak==0:
        return 0.,grad
    s=0.
    for i in range(n):
        for j in range(i):
            dot=0.
            for k in range(d):
                dot+=x[i,k]*x[j,k]
            v=(dot-target)/peak
            if v>0:
                w=v**(power-1)
                s+=w*v
                for k in range(d):
                    grad[i,k]+=w*x[j,k]
                    grad[j,k]+=w*x[i,k]
    grad*=2*peak*s**(2./power-1.)
    return peak*peak*s**(2./power),grad
