"""Shared utilities for the dimension-27 verification scripts.

Everything here is exact.  Coordinates live in the ring Z[sqrt2, sqrt3, sqrt6]
and are represented as integer 4-tuples in the basis (1, sqrt2, sqrt3, sqrt6).
Comparisons that cannot be settled by integer arithmetic alone are settled with
certified rational enclosures of the radicals, computed from integer square
roots.  No floating-point arithmetic is used anywhere in this module.

No third-party packages are required.
"""
import collections
import hashlib
import math
import re
from fractions import Fraction

# ---------------------------------------------------------------------------
# the ring Z[sqrt2, sqrt3, sqrt6]
# ---------------------------------------------------------------------------
ZERO = (0, 0, 0, 0)
_RAD_OF_SLOT = (1, 2, 3, 6)
_SLOT_OF_RAD = {1: 0, 2: 1, 3: 2, 6: 3}
_NAME = (None, "sqrt(2)", "sqrt(3)", "sqrt(6)")


def radd(x, y):
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])


def rsub(x, y):
    return (x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3])


def rmul(x, y):
    """Product in Z[sqrt2,sqrt3,sqrt6]: sqrt2*sqrt3 = sqrt6, sqrt2*sqrt6 = 2 sqrt3, ..."""
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return (a1 * a2 + 2 * b1 * b2 + 3 * c1 * c2 + 6 * d1 * d2,
            a1 * b2 + b1 * a2 + 3 * (c1 * d2 + d1 * c2),
            a1 * c2 + c1 * a2 + 2 * (b1 * d2 + d1 * b2),
            a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2)


def rdot(u, v):
    """Dot product of two sequences of ring elements."""
    acc = ZERO
    for a, b in zip(u, v):
        acc = radd(acc, rmul(a, b))
    return acc


def is_rational(x):
    return x[1] == 0 and x[2] == 0 and x[3] == 0


def rstr(x):
    """Render a ring element in the data set's token syntax."""
    a, b, c, d = (int(t) for t in x)
    if a == b == c == d == 0:
        return "0"
    out = str(a) if a else ""
    for slot, coef in ((1, b), (2, c), (3, d)):
        if not coef:
            continue
        sign = "-" if coef < 0 else ("+" if out else "")
        mag = abs(coef)
        out += sign + (f"{mag}*{_NAME[slot]}" if mag != 1 else _NAME[slot])
    return out


# re.ASCII matters: without it \d also matches non-ASCII decimal digits, so a
# coordinate written with, say, U+FF12 FULLWIDTH DIGIT TWO would parse as 2.
_TOKEN = re.compile(r"([+-]?)\s*(?:([0-9]+)\s*\*\s*)?sqrt\(([0-9]+)\)"
                    r"|([+-]?\s*[0-9]+)", re.ASCII)


def parse_token(tok):
    """'2*sqrt(3)+sqrt(6)' -> (0,0,2,1).

    Strict: every term after the first must carry an explicit sign, so
    juxtapositions such as '2sqrt(3)' or '2 3' are rejected rather than being
    silently read as products or sums.  Only ASCII digits are accepted.
    """
    s = tok.strip()
    if not s:
        raise ValueError("empty coordinate")
    if not s.isascii():
        raise ValueError(f"non-ASCII character in coordinate {tok!r}")
    pos = 0
    out = [0, 0, 0, 0]
    for m in _TOKEN.finditer(s):
        if m.start() != pos:
            raise ValueError(f"unparsed text in coordinate {tok!r}")
        lead = (m.group(1) or "") if m.group(4) is None else m.group(4).strip()[:1]
        # `lead not in "+-"` would be wrong: "" is a substring of every string,
        # so an unsigned term would slip through and '2sqrt(3)' would parse.
        if pos > 0 and lead not in ("+", "-"):
            raise ValueError(f"missing operator in coordinate {tok!r}")
        pos = m.end()
        if m.group(4) is not None:
            out[0] += int(m.group(4).replace(" ", ""))
        else:
            sign = -1 if m.group(1) == "-" else 1
            coef = int(m.group(2)) if m.group(2) else 1
            rad = int(m.group(3))
            if rad not in _SLOT_OF_RAD:
                raise ValueError(f"unsupported radical sqrt({rad}) in {tok!r}")
            out[_SLOT_OF_RAD[rad]] += sign * coef
    if pos != len(s):
        raise ValueError(f"unparsed text at end of coordinate {tok!r}")
    return tuple(out)


# ---------------------------------------------------------------------------
# certified comparison of ring elements
# ---------------------------------------------------------------------------
def _sqrt_bounds(n, den):
    r = math.isqrt(n * den * den)
    return Fraction(r, den), Fraction(r + 1, den)


def ring_bounds(x, den):
    """Rigorous [lo, hi] enclosing a + b sqrt2 + c sqrt3 + d sqrt6.

    The rational part x[0] may be an int or a Fraction; the rest are ints.
    """
    lo = hi = Fraction(x[0])
    for slot in (1, 2, 3):
        coef = int(x[slot])
        if not coef:
            continue
        a, b = _sqrt_bounds(_RAD_OF_SLOT[slot], den)
        lo += coef * (a if coef > 0 else b)
        hi += coef * (b if coef > 0 else a)
    return lo, hi


def ring_cmp_rational(x, q):
    """Sign of x - q for a ring element x and a Fraction q.  Exact.

    If x - q is rational the comparison is immediate.  Otherwise x - q is
    irrational, hence nonzero, and a tightening enclosure decides the sign.
    """
    diff = (Fraction(x[0]) - q, x[1], x[2], x[3])
    if diff[1] == diff[2] == diff[3] == 0:
        return (diff[0] > 0) - (diff[0] < 0)
    den = 10 ** 40
    while True:
        lo, hi = ring_bounds(diff, den)
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        den *= 10 ** 20          # x - q is irrational, so this terminates


def ring_cmp(x, y):
    """Sign of x - y for two ring elements.  Exact."""
    return ring_cmp_rational(rsub(x, y), Fraction(0))


def squarefree_part(n):
    """n = m^2 * d with d squarefree; returns (m, d).  n > 0."""
    m, d, k = 1, n, 2
    while k * k <= d:
        while d % (k * k) == 0:
            d //= k * k
            m *= k
        k += 1
    return m, d


def floor_half_sqrt_minus(n, x):
    """floor( sqrt(n)/2 - x ) exactly, for an integer n > 0 and ring element x.

    sqrt(n) lies in Q(sqrt2,sqrt3) exactly when the squarefree part of n is 1,
    2, 3 or 6; then 2*(sqrt(n)/2 - x) is a ring element and the value is either
    a half-integer or irrational.  Otherwise sqrt(n) is not in the field, so the
    value is irrational.  Irrational values are never integers, so a tightening
    enclosure determines the floor.
    """
    m, d = squarefree_part(n)
    twice = None
    if d in _SLOT_OF_RAD:
        twice = [-2 * int(t) for t in x]
        twice[_SLOT_OF_RAD[d]] += m
        if twice[1] == twice[2] == twice[3] == 0:
            return twice[0] // 2                     # exactly a half-integer
    den = 10 ** 40
    while True:
        if twice is not None:
            lo, hi = ring_bounds(tuple(twice), den)
            lo, hi = lo / 2, hi / 2
        else:
            a, b = _sqrt_bounds(n, den)
            xlo, xhi = ring_bounds(x, den)
            lo, hi = a / 2 - xhi, b / 2 - xlo
        f_lo, f_hi = math.floor(lo), math.floor(hi)
        if f_lo == f_hi:
            return f_lo
        den *= 10 ** 20


# ---------------------------------------------------------------------------
# GF(2) linear algebra on 24-bit words
# ---------------------------------------------------------------------------
def gf2_span(generators):
    """Return (basis, set of all codewords) of the span of the given words.

    Gaussian elimination keyed by leading bit: each basis vector owns a distinct
    leading bit, so the span has exactly 2**len(basis) elements.
    """
    pivot = {}
    for g in generators:
        cur = g
        while cur:
            top = cur.bit_length() - 1
            if top not in pivot:
                pivot[top] = cur
                break
            cur ^= pivot[top]
    basis = [pivot[k] for k in sorted(pivot, reverse=True)]
    code = {0}
    for b in basis:
        code |= {c ^ b for c in code}
    return basis, code


def weight_distribution(code):
    return dict(sorted(collections.Counter(bin(c).count("1") for c in code).items()))


def support_mask(vec):
    m = 0
    for i, v in enumerate(vec):
        if v:
            m |= 1 << i
    return m


# ---------------------------------------------------------------------------
# the Leech lattice
# ---------------------------------------------------------------------------
def in_leech(x, golay):
    """Is the integer vector x in Lambda_24, in the scaling where the minimal
    vectors have squared norm 32?

    Conway and Sloane, SPLAG chapter 10: x lies in the lattice exactly when all
    coordinates are congruent to a common m modulo 2, the set of positions with
    x_i congruent to m + 2 modulo 4 is a Golay codeword, and the coordinate sum
    is congruent to 4m modulo 8.
    """
    m = x[0] & 1
    s = 0
    total = 0
    for i, v in enumerate(x):
        if (v & 1) != m:
            return False
        total += v
        if (v - m - 2) % 4 == 0:
            s |= 1 << i
    if s not in golay:
        return False
    return total % 8 == (4 * m) % 8


# the three shapes of minimal vector of Lambda_24, written as the multiset of
# absolute values of the nonzero coordinates: (+-4^2, 0^22), (+-2^8, 0^16) and
# (-+3, +-1^23).
LEECH_SHAPE_CENSUS = {((4, 2),): 1104,
                      ((2, 8),): 97152,
                      ((3, 1), (1, 23)): 98304}
GOLAY_WEIGHT_DISTRIBUTION = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}


def leech_shape(vec):
    """The multiset of absolute values of the nonzero coordinates,
    as a tuple of (value, multiplicity) pairs in decreasing order of value."""
    counts = {}
    for v in vec:
        if v:
            a = abs(v)
            counts[a] = counts.get(a, 0) + 1
    return tuple(sorted(counts.items(), reverse=True))


def certify_leech_minimal(heads):
    """Facts certifying that `heads` are minimal vectors of a copy of Lambda_24.

    Returns a dict; the caller turns the entries into checks.  Nothing here is
    assumed about the input beyond its being a list of integer 24-vectors.
    """
    census = {}
    for y in heads:
        s = leech_shape(y)
        census[s] = census.get(s, 0) + 1
    octads = sorted({support_mask(y) for y in heads if leech_shape(y) == ((2, 8),)})
    basis, golay = gf2_span(octads)
    bad = next((y for y in heads if not in_leech(y, golay)), None)
    return {
        "count": len(heads),
        "distinct": len(set(heads)) == len(heads),
        "norms_32": all(sum(x * x for x in y) == 32 for y in heads),
        "census": census,
        "census_ok": census == LEECH_SHAPE_CENSUS,
        "n_octads": len(octads),
        "golay_dim": len(basis),
        "golay_size": len(golay),
        "weights": weight_distribution(golay),
        "weights_ok": weight_distribution(golay) == GOLAY_WEIGHT_DISTRIBUTION,
        "octads_are_weight8":
            set(octads) == {c for c in golay if bin(c).count("1") == 8},
        "all_in_leech": bad is None,
        "first_bad": bad,
    }


# ---------------------------------------------------------------------------
# splitting a dimension-27 configuration into its three families
# ---------------------------------------------------------------------------
def canonical_scale(norm, target):
    """(p, q) meaning "multiply the row by p/q" to send squared norm to target.

    p/q = sqrt(target*norm)/norm, which is rational exactly when target*norm is
    a perfect square.  Returns None otherwise.
    """
    prod = target * norm
    root = math.isqrt(prod)
    return (root, norm) if root * root == prod else None


def _rescale_int(vec, num, den, what):
    out = []
    for v in vec:
        v *= num
        if v % den:
            raise ValueError(f"{what} does not rescale to an integer vector")
        out.append(v // den)
    return tuple(out)


def _rescale_ring(vec, num, den, what):
    """The same, for a sequence of ring elements."""
    return tuple(_rescale_int(coord, num, den, what) for coord in vec)


def _scale_to(total, target, what):
    sc = canonical_scale(total, target)
    if sc is None:
        raise ValueError(f"{what}: cannot be rescaled to squared norm {target}")
    return sc


def split_families(cfg):
    """Split a dimension-27 configuration into equator, cap and axis families.

    Each row is rescaled to a canonical length: equator rows to squared norm 32,
    cap rows to 432 (head 288 = 9*32, tail 144) and axis rows to 48.  Every step
    is guarded -- the rescaling must be by a rational whose square is exact, the
    result must be integral, and a cap head must be three times an integer
    vector.  Anything else raises ValueError rather than being waved through.

    Returns (equator, caps, axis, rational_columns, radical_columns) where
    equator is a list of integer 24-tuples y with y.y = 32, caps a list of
    (y, t) with y.y = 32 and t.t = 144, and axis a list of t with t.t = 48.
    """
    dim = cfg.dimension
    if dim != 27:
        raise ValueError(f"expected dimension 27, got {dim}")
    rational_col = [all(is_rational(row[j]) for row in cfg.rows) for j in range(dim)]
    rat = [j for j in range(dim) if rational_col[j]]
    irr = [j for j in range(dim) if not rational_col[j]]

    equator, caps, axis = [], [], []
    for r, row in enumerate(cfg.rows):
        head = tuple(row[j][0] for j in rat)
        tail = tuple(row[j] for j in irr)
        hn = sum(x * x for x in head)
        tn = rdot(tail, tail)
        if not is_rational(tn):
            raise ValueError(f"row {r}: irrational squared norm")
        tn = tn[0]
        total = hn + tn
        if total <= 0:
            raise ValueError(f"row {r}: zero row")
        where = f"row {r}"
        if tn == 0:                                   # equator
            num, den = _scale_to(total, 32, where)
            equator.append(_rescale_int(head, num, den, where))
        elif hn == 0:                                 # axis
            num, den = _scale_to(total, 48, where)
            axis.append(_rescale_ring(tail, num, den, where))
        elif 3 * tn == total:                         # cap
            num, den = _scale_to(total, 432, where)
            h = _rescale_int(head, num, den, where)
            t = _rescale_ring(tail, num, den, where)
            if any(v % 3 for v in h):
                raise ValueError(f"{where}: cap head is not 3 times an integer vector")
            caps.append((tuple(v // 3 for v in h), t))
        else:
            raise ValueError(f"{where}: |w|^2/|p|^2 is neither 0, 1/3 nor 1")
    return equator, caps, axis, rat, irr


def shape_str(shape):
    return " ".join(f"(+-{v})^{m}" for v, m in shape)


# ---------------------------------------------------------------------------
# configuration files
# ---------------------------------------------------------------------------
def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class Configuration:
    """One 'Dimension: n' block of a kissing-number data file, parsed exactly."""

    def __init__(self, dimension, coefficients, rows):
        self.dimension = dimension
        self.coefficients = coefficients
        self.rows = rows                     # list of tuples of ring elements

    def __len__(self):
        return len(self.rows)


def parse_configuration(path, expect_dimension=None):
    with open(path) as f:
        lines = f.read().split("\n")
    while lines and lines[-1] == "":
        lines.pop()
    if len(lines) < 5:
        raise ValueError("file is too short to be a configuration block")
    if not lines[0].startswith("Dimension:"):
        raise ValueError("first line must be 'Dimension: n'")
    dimension = int(lines[0].split(":", 1)[1])
    if not lines[1].startswith("Number of points:"):
        raise ValueError("second line must be 'Number of points: N'")
    count = int(lines[1].split(":", 1)[1])
    if not lines[2].startswith("Inner product coefficients:"):
        raise ValueError("third line must be 'Inner product coefficients: ...'")
    coefficients = [int(t) for t in lines[2].split(":", 1)[1].strip().split(",")]
    if lines[3].strip() != "Points:":
        raise ValueError("fourth line must be 'Points:'")
    body = lines[4:]
    if len(body) != count:
        raise ValueError(f"header declares {count} points, file has {len(body)} rows")
    if len(coefficients) != dimension:
        raise ValueError("wrong number of inner product coefficients")
    if expect_dimension is not None and dimension != expect_dimension:
        raise ValueError(f"expected dimension {expect_dimension}, file says {dimension}")

    cache = {}
    rows = []
    for r, line in enumerate(body):
        parts = line.split(",")
        if len(parts) != dimension:
            raise ValueError(f"row {r} has {len(parts)} coordinates, expected {dimension}")
        row = []
        for p in parts:
            v = cache.get(p)
            if v is None:
                v = cache[p] = parse_token(p)
            row.append(v)
        rows.append(tuple(row))
    cfg = Configuration(dimension, coefficients, rows)
    cfg.tokens = set(cache)
    return cfg


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
class Report:
    """Collects pass/fail lines and prints them as they happen."""

    def __init__(self, title):
        self.failures = 0
        self.checks = 0
        print("=" * 78)
        print(title)
        print("=" * 78, flush=True)

    def check(self, ok, label, detail=""):
        self.checks += 1
        if not ok:
            self.failures += 1
        mark = "ok  " if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f"   {detail}" if detail else ""), flush=True)
        return ok

    def note(self, text):
        print(f"         {text}", flush=True)

    def section(self, text):
        print(f"\n{text}", flush=True)

    def done(self, success_message):
        print()
        if self.failures == 0:
            print(f"{self.checks}/{self.checks} checks passed.")
            print(success_message, flush=True)
            return 0
        print(f"{self.checks - self.failures}/{self.checks} checks passed, "
              f"{self.failures} FAILED.", flush=True)
        return 1
