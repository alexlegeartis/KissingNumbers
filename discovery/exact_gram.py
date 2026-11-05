"""Independent bounded rational Gram certification and exact integer span rank.

Rational JSON scalars are integers or [numerator, positive denominator].
No floating point values or serialized proof flags enter these checks.
"""
from fractions import Fraction
from math import gcd
from pathlib import Path

from .core import RunContext, Verification, digest


def bounded_fraction(value, context, label='rational scalar'):
    if isinstance(value, bool):
        raise ValueError(f'{label} must be an integer or [numerator, denominator]')
    if isinstance(value, int):
        context.limit('max_rational_bits', value.bit_length())
        return Fraction(value)
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f'{label} must be an integer or [numerator, denominator]')
    numerator, denominator = value
    if any(isinstance(x, bool) or not isinstance(x, int) for x in value) or denominator <= 0:
        raise ValueError(f'{label} requires integer numerator and positive integer denominator')
    # Bound raw integers before Fraction computes a gcd or multiplies anything.
    context.limit('max_rational_bits', numerator.bit_length())
    context.limit('max_rational_bits', denominator.bit_length())
    return Fraction(numerator, denominator)


def rational_json(value):
    value = Fraction(value)
    return value.numerator if value.denominator == 1 else [value.numerator, value.denominator]


def check_fraction_size(value, context):
    context.limit('max_rational_bits', max(value.numerator.bit_length(), value.denominator.bit_length()))
    return value


def exact_integer_rank(rows, context):
    """Fraction-free row elimination; count scalar updates before performing them.

    Dividing a row by its gcd preserves its rational span and controls growth.
    At most one pivot row per coordinate is retained.
    """
    if not rows:
        return 0, 0
    n, width = len(rows), len(rows[0])
    context.limit('max_materialized_points', n)
    context.limit('max_dimension', width)
    # Validate every row even when an early full-rank witness exists.
    for source in rows:
        if not isinstance(source, list) or len(source) != width:
            raise ValueError('Ragged integer rank matrix')
        for value in source:
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError('Exact rank input must contain integers')
            context.limit('max_rational_bits', value.bit_length())
    basis = {}; work = 0
    for source in rows:
        row = source[:]
        for pivot, old in basis.items():
            if not row[pivot]:
                continue
            work += 2*width
            context.limit('max_rank_work', work)
            coefficient, old_pivot = row[pivot], old[pivot]
            row = [old_pivot*x-coefficient*y for x, y in zip(row, old)]
            common = 0
            for value in row:
                context.limit('max_rational_bits', value.bit_length())
                common = gcd(common, value)
            if common > 1:
                row = [value//common for value in row]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            common = 0
            for value in row[pivot:]:
                common = gcd(common, value)
            if row[pivot] < 0:
                common = -common
            basis[pivot] = [value//common for value in row]
            # Input rows need not arrive with increasing pivots. Preserve a
            # sorted echelon basis so every new row is reduced in pivot order.
            basis = dict(sorted(basis.items()))
            # A full-width independent set is an exact upper and lower bound;
            # later rows cannot change rank and have already been validated.
            if len(basis)==width:return width,work
    return len(basis), work


def positive_semidefinite_rank(matrix, context):
    """Exact symmetric Schur elimination, including singular residual checks."""
    residual = [row[:] for row in matrix]
    rank = 0; work = 0
    while residual:
        n = len(residual)
        negative = next((i for i in range(n) if residual[i][i] < 0), None)
        if negative is not None:
            raise ValueError('Gram matrix is not positive semidefinite (negative Schur diagonal)')
        pivot = next((i for i in range(n) if residual[i][i] > 0), None)
        if pivot is None:
            if any(value != 0 for row in residual for value in row):
                raise ValueError('Gram matrix is indefinite (zero diagonal with nonzero residual)')
            break
        if pivot:
            residual[0], residual[pivot] = residual[pivot], residual[0]
            for row in residual:
                row[0], row[pivot] = row[pivot], row[0]
        context.limit('max_dimension', rank+1)
        diagonal = residual[0][0]
        tail = n-1
        # Each off-diagonal update uses division, multiplication, subtraction;
        # count all three, including the diagonal updates, before allocation.
        work += 3*tail*(tail+1)//2
        context.limit('max_rank_work', work)
        ratios = [check_fraction_size(residual[i][0]/diagonal, context) for i in range(1, n)]
        reduced = [[Fraction(0) for _ in range(tail)] for _ in range(tail)]
        for i in range(tail):
            for j in range(i, tail):
                value = check_fraction_size(residual[i+1][j+1]-ratios[i]*residual[0][j+1], context)
                reduced[i][j] = reduced[j][i] = value
        residual = reduced
        rank += 1
        context.limit('max_dimension', rank)
    return rank, work


def verify_rational_gram(data, context=None):
    context = context or RunContext(Path('.'))
    if not isinstance(data, dict) or set(data)-{'gram', 'norm_squared'}:
        raise ValueError('Gram data requires gram and optional norm_squared only')
    raw = data.get('gram')
    if not isinstance(raw, list) or not raw:
        raise ValueError('Gram matrix must be a nonempty square list')
    n = len(raw)
    context.limit('max_materialized_points', n)
    context.limit('max_gram_entries', n*n)
    context.limit('max_pairs', n*(n-1)//2)
    if any(not isinstance(row, list) or len(row) != n for row in raw):
        raise ValueError('Gram matrix must be square')
    matrix = [[bounded_fraction(value, context, 'Gram entry') for value in row] for row in raw]
    norm = matrix[0][0]
    if norm <= 0:
        raise ValueError('Common Gram diagonal must be positive')
    if 'norm_squared' in data and bounded_fraction(data['norm_squared'], context, 'norm_squared') != norm:
        raise ValueError('Claimed norm does not match Gram diagonal')
    maximum = None
    for i in range(n):
        if matrix[i][i] != norm:
            raise ValueError('Gram diagonal is not common')
        for j in range(i):
            value = matrix[i][j]
            if value != matrix[j][i]:
                raise ValueError('Gram matrix is not symmetric')
            if value == norm:
                raise ValueError('Duplicate points in Gram matrix')
            if 2*value > norm:
                raise ValueError('Gram pairwise kissing inequality violated')
            maximum = value if maximum is None else max(maximum, value)
    rank, work = positive_semidefinite_rank(matrix, context)
    facts = {'dimension': rank, 'rank_exact': rank, 'ambient_dimension': rank,
             'dimension_scope': 'exact Gram rank; symbolic real realization',
             'cardinality': n, 'norm_squared': rational_json(norm),
             'max_inner_product': None if maximum is None else rational_json(maximum),
             'pairs_checked': n*(n-1)//2}
    return Verification('exact', 'exact rational symmetric Gram, Schur PSD/rank, common norm and all pair inequalities',
                        {'facts': facts, **facts, 'rank_scalar_work': work, 'data_sha256': digest(data),
                         'realization': 'A PSD rational Gram matrix has a real realization in its exact rank; no rational coordinates are asserted.'})
