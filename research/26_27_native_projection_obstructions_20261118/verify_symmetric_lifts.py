"""Reproduce exact linear-programming bounds for symmetric one-axis lifts."""

from fractions import Fraction as F
import json


def gegenbauer(dimension, t, degree):
    """Normalized Gegenbauer polynomials evaluated at a rational number."""
    q = [F(1), t]
    for k in range(1, degree):
        a = dimension - 2
        q.append(((2 * k + a) * t * q[k] - k * q[k - 1]) / (k + a))
    return q


def evaluate_section24():
    # These are the distinct products <= 1/3 in the normalized Leech shell.
    allowed = [F(-1), F(-1, 2), F(-1, 4), F(0), F(1, 4)]
    degrees = [1, 2, 3, 4, 9]
    coefficients = [
        F(928, 77),
        F(5083, 77),
        F(25826056, 129591),
        F(27600, 77),
        F(27747200, 129591),
    ]
    if min(coefficients) < 0:
        raise ValueError("Negative Gegenbauer coefficient")
    for t in allowed:
        q = gegenbauer(24, t, 9)
        if 1 + sum(c * q[i] for i, c in zip(degrees, coefficients)) > 0:
            raise ValueError(f"Positive value at {t}")
    return 1 + sum(coefficients)


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def evaluate_sphere26():
    roots = [F(-9, 14), F(-5, 18), F(1, 38)]
    h = [F(1)]
    for x in roots:
        h = multiply(h, [-x, F(1)])
    # f(t) = (t - 1/3) h(t)^2 has the required sign for every t <= 1/3.
    f = multiply([-F(1, 3), F(1)], multiply(h, h))
    q = [[F(1)], [F(0), F(1)]]
    a = 26 - 2
    for k in range(1, 7):
        # Keep this arithmetic in Q, including the division by k + a.
        nxt = [F(0)] + [F(2 * k + a, k + a) * x for x in q[k]]
        for j, x in enumerate(q[k - 1]):
            nxt[j] -= F(k, k + a) * x
        q.append(nxt)
    remainder = f.copy()
    coeff = [F(0)] * len(f)
    for i in reversed(range(len(f))):
        coeff[i] = remainder[i] / q[i][-1]
        for j, x in enumerate(q[i]):
            remainder[j] -= coeff[i] * x
    if any(remainder) or coeff[0] <= 0 or min(coeff) < 0:
        raise ValueError("Invalid Gegenbauer expansion")
    return sum(f) / coeff[0]


def main():
    leech_upper = evaluate_section24()
    sphere_upper = evaluate_sphere26()
    # At most one extra split uses the two poles.
    result = {
        "normalized_Leech_subset_LP": str(leech_upper),
        "sphere_26_subset_LP": str(sphere_upper),
        "symmetric_25_to_26_upper": 197586 + int(leech_upper) + 2 * 926 + 247 + 1,
        "symmetric_26_to_27_upper": 201192 + int(sphere_upper) + 1,
    }
    if result["symmetric_25_to_26_upper"] > 201192:
        raise ValueError("The lift bound misses the 26D target")
    if result["symmetric_26_to_27_upper"] > 210108:
        raise ValueError("The lift bound misses the 27D target")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
