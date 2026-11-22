"""Exact complete-fibre calculation for a norm-eight normal in an extremal P48."""
from fractions import Fraction
from itertools import combinations
from math import prod

import sympy as sp


def layer_sizes(norm, shell_size, max_inner, exceptional):
    # For the norm-eight shell remove the two exceptional vectors ±w first.
    moments = []
    for j in range(max_inner+1):
        mean = Fraction((8*norm)**j * prod(range(1, 2*j, 2)), prod(range(48, 48+2*j, 2)))
        moments.append(shell_size*mean - exceptional*8**(2*j))
    vandermonde = sp.Matrix([[i**(2*j) for i in range(max_inner+1)] for j in range(max_inner+1)])
    both = vandermonde.inv() * sp.Matrix(moments)
    assert all(x.is_Integer and x >= 0 for x in both)
    assert all(int(x) % 2 == 0 for x in both[1:])
    return [int(both[0]), *(int(x)//2 for x in both[1:])]


def main():
    minimum = layer_sizes(6, 52_416_000, 4, 0)
    # The coefficient of q^4 counts the norm-eight shell.
    q = sp.Symbol("q")
    e4 = sum(240*sp.divisor_sigma(i, 3)*q**i for i in range(1, 5)) + 1
    delta = q-24*q**2+252*q**3-1472*q**4
    expansion = sp.series(e4**6 - 1440*e4**3*delta + 125280*delta**2, q, 0, 5).removeO()
    assert expansion.coeff(q,1) == expansion.coeff(q,2) == 0 and expansion.coeff(q,3) == 52_416_000
    next_shell = layer_sizes(8, int(expansion.coeff(q, 4)), 5, 2)
    # In each congruence class of (x,w) modulo eight choose the representative
    # with |(x,w)| <= 4. Projected squared radii are in units of 1/8.
    families = [(j % 8, 48-j*j, minimum[abs(j)]) for j in range(-3, 5)]
    families.append((4, 48, next_shell[4]))
    distance = [48, 47, 44, 39, 32, 39, 44, 47]
    best = (0, ())
    for k in range(1, len(families)+1):
        for config in combinations(families, k):
            if all((p := r+s-distance[(i-j)%8]) <= 0 or p*p <= r*s
                   for (i,r,_),(j,s,_) in combinations(config, 2)):
                best = max(best, (sum(w for _,_,w in config), config))
    assert best[0] == 20_876_648
    # With one real head and total squared radius <= 6, every family except
    # the radius-four family permits one copy: m(2r-6) <= 6. The latter
    # permits two copies on a line (its latitude angle is < 45 degrees).
    other = sum(n for _,r,n in families if r != 48)
    small = minimum[4]
    assert minimum[0] + other + small == 52_416_000
    assert next_shell[4] + other + small + 2 < 52_416_000
    assert other + small + 2 < 52_416_000
    assert minimum[0] + other + 2 < 52_416_000
    print("47D maximum under the class-distance test:", best)
    print("48D maximum under the one-dimensional head and class-distance constraints: 52416000")


if __name__ == "__main__":
    main()
