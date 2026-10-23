"""Exact rational lattice enumeration; standard library only."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import time

HERE = Path(__file__).resolve().parent


def transpose(A):
    return list(map(list, zip(*A)))


def multiply(A, B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def determinant(A):
    A = [[F(x) for x in row] for row in A]
    value = F(1)
    for j in range(len(A)):
        pivot = next((i for i in range(j,len(A)) if A[i][j]),None)
        if pivot is None:
            return F(0)
        if pivot != j:
            A[j], A[pivot] = A[pivot], A[j]
            value = -value
        v = A[j][j]
        value *= v
        for i in range(j+1,len(A)):
            factor = A[i][j]/v
            for k in range(j+1,len(A)):
                A[i][k] -= factor*A[j][k]
    return value


def enumerate_short(G, bound=8):
    n = len(G)
    L = [[F(i==j) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        D.append(F(G[j][j])-sum(L[j][k]**2*D[k] for k in range(j)))
        assert D[-1] > 0, 'Gram matrix is not positive definite'
        for i in range(j+1,n):
            L[i][j] = (G[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/D[j]
    v = [0]*n
    def visit(k, remaining):
        if k < 0:
            norm = bound-remaining
            assert norm.denominator == 1
            if norm:
                yield tuple(v), int(norm)
            return
        shift = sum((L[i][k]*v[i] for i in range(k+1,n)),F(0))
        r2 = remaining/D[k]
        # Strictly larger than sqrt(r2), computed using integer arithmetic.
        radius = math.isqrt(r2.numerator//r2.denominator)+1
        center = (-shift.numerator)//shift.denominator
        for z in range(center-radius-1, center+radius+2):
            cost = D[k]*(F(z)+shift)**2
            if cost <= remaining:
                v[k] = z
                yield from visit(k-1, remaining-cost)
    yield from visit(n-1,F(bound))


