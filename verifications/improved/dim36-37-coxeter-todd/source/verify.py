"""Standalone exact verifier for a 36-dimensional Coxeter--Todd coset code.

Python 3.9+, standard library only. Run from any directory:
    python path/to/verify.py
The verifier enumerates all K12 vectors of norm <= 8 using rational LDL,
then checks every pair of the explicitly supplied outer cosets.
It does not trust the search, its quadratic form, its LPs, or its solvers.
"""
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


def main():
    started = time.monotonic()
    path = HERE/'certificate.json'
    cert = json.loads(path.read_text())
    assert cert['dimension'] == 36
    assert cert['squared_norm'] == 12
    G, M = cert['gram_12'], cert['sublattice_matrix']
    assert len(G) == 12 and all(len(row)==12 for row in G)
    assert G == transpose(G)
    assert all(isinstance(x,int) for row in G+M for x in row)
    assert multiply(multiply(M,G),transpose(M)) == [[3*x for x in row] for row in G]
    assert abs(determinant(M)) == 729
    first, second = cert['label_first'], cert['label_second']
    assert len(first) == len(second) == 6
    assert sorted(first+second) == list(range(12))
    powers = [3**i for i in range(6)]
    def label(v):
        return sum(((v[i]+v[j])%3)*p for i,j,p in zip(first,second,powers))
    assert all(label(row)==0 for row in M)
    # The label map is onto: the six unit vectors at `first` give its basis.
    # Its kernel therefore has index 729, equal to the row lattice of M.
    short = list(enumerate_short(G))
    shell_counts = Counter(norm for _,norm in short)
    assert shell_counts == {4:756,6:4032,8:20412}, shell_counts
    minima = [100]*729
    multiplicities = [0]*729
    minima[0], multiplicities[0] = 0,1
    for v,norm in short:
        assert sum(v[i]*G[i][j]*v[j] for i in range(12) for j in range(12)) == norm
        c = label(v)
        if norm < minima[c]:
            minima[c], multiplicities[c] = norm,1
        elif norm == minima[c]:
            multiplicities[c] += 1
    assert Counter(zip(minima,multiplicities)) == {(0,1):1,(4,3):252,(6,18):224,(8,81):252}
    # Exhaustive enumeration proves these coset minima, with no polynomial
    # classification or floating calculation assumed. M K12 has minimum 12
    # and 756 minimal vectors by the verified similarity M G M^T=3G.
    outer = cert['outer_cosets']
    assert outer and outer[0] == [0,0,0]
    assert all(len(row)==3 and all(isinstance(x,int) and 0<=x<729 for x in row) for row in outer)
    assert len(set(map(tuple,outer))) == len(outer)
    assert all(sum(minima[x] for x in row)==12 for row in outer[1:])
    total = 3*shell_counts[4] + sum(math.prod(multiplicities[x] for x in row) for row in outer[1:])
    assert total == cert['points'], (total,cert['points'])
    digits = [[(i//p)%3 for p in powers] for i in range(729)]
    distance = [[minima[sum(((a-b)%3)*p for a,b,p in zip(x,y,powers))]
                 for y in digits] for x in digits]
    smallest = 100
    checked = 0
    for i,(a,b,c) in enumerate(outer):
        da,db,dc = distance[a],distance[b],distance[c]
        for j in range(i):
            x,y,z = outer[j]
            gap = da[x]+db[y]+dc[z]
            assert gap >= 12, ('conflicting cosets',i,j,gap)
            if gap < smallest:
                smallest = gap
            checked += 1
    # Distinct points in the same coset differ by a nonzero vector of
    # (M K12)^3, also of norm at least 12. All points have norm exactly 12.
    result = {'dimension':36,'points':total,'squared_norm':12,
              'minimum_squared_distance_at_least':12,'maximum_inner_product_at_most':6,
              'unit_sphere_inner_product_at_most':'1/2','outer_cosets':len(outer),
              'outer_pairs_checked':checked,'minimum_outer_distance':smallest,
              'exact_short_shell_counts':dict(shell_counts),'all_checks_pass':True,
              'arithmetic':'Python integers and fractions.Fraction; no third-party libraries',
              'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
              'elapsed_seconds':time.monotonic()-started}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
