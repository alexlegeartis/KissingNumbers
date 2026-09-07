#!/usr/bin/env python3
"""Fact-check every arithmetic and cross-reference claim in kissing46.tex against the
repository's own single sources of truth.  Recomputes rather than parses: the moment
identities, the cap count, Caro-Wei, the Edel-Rains-Sloane totals, the class bound and
every entry of Tables 1, 3 and 4 are derived here and compared with what the paper says.

    python factcheck.py        # a second, no third-party packages
"""
import io
import json, math, os, re, sys
from fractions import Fraction as F
from math import comb, floor, ceil, sqrt

HERE = os.path.dirname(os.path.abspath(__file__))


def _find_repo():
    """The verification repository, whether or not this file sits inside it."""
    here = os.path.dirname(HERE)
    for cand in (here, os.path.join(here, 'kissing_verifications')):
        if os.path.exists(os.path.join(cand, 'RESULTS.md')):
            return cand
    raise SystemExit('cannot find the verification repository: expected RESULTS.md in\n'
                     '  %s\nor\n  %s' % (here, os.path.join(here, 'kissing_verifications')))


KV = _find_repo()
sys.path.insert(0, os.path.join(KV, 'common'))
from published import COHN, TAU, floor_for

TEX = io.open(os.path.join(HERE, 'kissing46.tex'), encoding='utf-8').read()

def _npy_rows(path):
    """Rows of a .npy, from its header alone -- factcheck.py has no numpy."""
    import ast
    with io.open(path, 'rb') as h:
        head = h.read(256)
    i = head.index(b'{')
    j = head.index(b'}', i) + 1
    return ast.literal_eval(head[i:j].decode('latin-1'))['shape'][0]


def _npz_rows(path, key):
    """Rows of one member of a .npz, from its header alone -- still no numpy."""
    import ast
    import zipfile
    with zipfile.ZipFile(path) as z:
        with z.open(key + '.npy') as h:
            head = h.read(256)
    i = head.index(b'{')
    j = head.index(b'}', i) + 1
    return ast.literal_eval(head[i:j].decode('latin-1'))['shape'][0]


def _layer_fractions():
    """{k: (|Z|, poles)} for every dimension that ships a rotated axis layer."""
    import json
    import zipfile
    pkg = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps')
    data = os.path.join(pkg, 'data')
    rows = {r['k']: r for r in json.load(open(os.path.join(pkg, 'rows81-95.json')))}
    out = {}
    for k in range(3, 24):
        # The layer file follows the CAP SET: poles_kap_<k>.npz is a rotated copy of the
        # Kappa section at k = 12 and 13, and poles_rot_<k>.npz a copy of Lambda_k or of the
        # record configuration elsewhere.  Both may never be read at the same k, or |Z| and
        # the layer would be measured against different configurations.
        rot = [f for f in ('poles_kap_%d.npz' % k, 'poles_rot_%d.npz' % k,
                           'poles_rot_k%d.npz' % k)
               if os.path.exists(os.path.join(data, f))]
        if not rot:
            continue
        with zipfile.ZipFile(os.path.join(data, rot[0])) as z:
            with z.open('keep.npy') as h:
                head = h.read(128)
        i = head.index(b'(') + 1
        poles = int(head[i:head.index(b',', i)])
        if k <= 8:
            nz = TAU[k]              # the direction set is the root system, and D5, E6, E7,
        else:                        # E8 all attain tau(k); the paper says so in S5.6
            src = rows[k]['W']
            f = ('kap%d_W.npy' % k if src == 'kappa' else
                 'rec_%d_W.npy' % k if src.startswith('record') else
                 'lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k)
            nz = _npy_rows(os.path.join(data, f))
        out[k] = (nz, poles)
    return out



def _shipped_poles(k):
    """How many poles dimension 72+k actually ships, read from the package's data.

    Five files can hold a layer -- from the record configuration of R^k, from the same idea
    over Z[sqrt2], from a rotated copy of the cap directions under either naming, and since
    2026-08-30 from a rotated copy of the KAPPA cap directions at k = 12 and 13 -- and the
    dimension ships the LARGEST, which is how final81-95.py chooses.  Reading only the
    rotated file would check the paper against a layer the paper does not claim.

    poles_rot_<k>.npz and poles_kap_<k>.npz are never both present at one k, because they
    are copies of DIFFERENT cap sets and a dimension uses one cap set; the Lambda_k layers
    at k = 12 and 13 were retired when the Kappa sections replaced them.
    """
    import zipfile
    d = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
    got = []
    for f in ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k, 'poles_kap_%d.npz' % k,
              'poles_rot_%d.npz' % k, 'poles_rot_k%d.npz' % k):
        p = os.path.join(d, f)
        if os.path.exists(p):
            with zipfile.ZipFile(p) as z:          # no numpy: read the .npy header
                with z.open('keep.npy') as h:
                    head = h.read(128)
            i = head.index(b'(') + 1
            got.append(int(head[i:head.index(b',', i)]))
    if not got:
        raise SystemExit('no pole layer for k=%d' % k)
    return max(got)


_BS = chr(92) + chr(92)          # a literal backslash, for the LaTeX-macro regexes
FLAT = re.sub(r'\s+', ' ', TEX)
RES = io.open(os.path.join(KV, 'RESULTS.md'), encoding='utf-8').read()

BS = chr(92)
BAD = []
OK = 0


def chk(name, cond, detail=''):
    global OK
    if cond:
        OK += 1
    else:
        BAD.append((name, detail))
        print('  FAIL  %-52s %s' % (name, detail))


def num(s):
    return int(re.sub(r'[^\d]', '', s))


def _grp(v):
    """An integer as the paper prints it: thin spaces every three digits."""
    s, out = str(v), ''
    while len(s) > 3:
        out, s = BS + ',' + s[-3:] + out, s[:-3]
    return s + out


def _npy_ints(path):
    """A one-dimensional integer .npy, read with the standard library alone.

    factcheck.py advertises itself as needing no third-party package, and importing numpy
    for one array sum would quietly cost that.  The format is documented and trivial: a
    magic string, a version, a header length, a dict literal, then the raw buffer.
    """
    import ast, array
    with open(path, 'rb') as fh:
        assert fh.read(6) == b'\x93NUMPY', path
        major = ord(fh.read(1)); fh.read(1)
        n = int.from_bytes(fh.read(2 if major == 1 else 4), 'little')
        head = ast.literal_eval(fh.read(n).decode('latin1').strip())
        assert not head['fortran_order'] and len(head['shape']) == 1, head
        # class_sizes_gpu.npy is int16: the widths actually shipped, not just the
        # two that happened to be needed when this was written.
        code = {'|i1': 'b', '|u1': 'B', '<i2': 'h', '<u2': 'H',
                '<i4': 'i', '<u4': 'I', '<i8': 'q', '<u8': 'Q'}[head['descr']]
        a = array.array(code)
        a.frombytes(fh.read())
    assert len(a) == head['shape'][0], (len(a), head['shape'])
    return list(a)


# ---------------------------------------------------------------- Table 4 vs RESULTS.md
print('== Table 4 (longtable) against RESULTS.md and published.py ==')
tab4 = {}
body = TEX[TEX.index('\\endfoot'):TEX.index('\\end{longtable}')]
for line in body.split('\n'):
    cells = line.split('&')
    if len(cells) != 4 or not cells[0].strip().startswith('$'):
        continue
    try:
        d = int(re.sub(r'[^0-9]', '', cells[0]))
    except ValueError:
        continue
    tab4[d] = (num(cells[1]), num(cells[2]), cells[3].replace(chr(92) * 2, '').strip())

results = {}
for m in re.finditer(r'^\| (\d+) \| ([\d\u202f ]+) \| \*\*([\d\u202f ]+)\*\* \|', RES, re.M):
    results[int(m.group(1))] = (num(m.group(2)), num(m.group(3)))

chk('Table 4 row count == 47', len(tab4) == 47, 'got %d' % len(tab4))
chk('RESULTS.md row count == 47', len(results) == 47, 'got %d' % len(results))
chk('dimensions 46 and 47 are not counted', not ({46, 47} & set(tab4)),
    'still present: %s' % sorted({46, 47} & set(tab4)))
chk('but the paper still records them', '23\\,766\\,960' in TEX
    and '12\\,309\\,600' in TEX)
chk('Table 4 dims == RESULTS dims', set(tab4) == set(results),
    'tex only %s / results only %s' % (sorted(set(tab4) - set(results)),
                                       sorted(set(results) - set(tab4))))
for d in sorted(set(tab4) & set(results)):
    chk('dim %d prev matches RESULTS.md' % d, tab4[d][0] == results[d][0],
        '%d vs %d' % (tab4[d][0], results[d][0]))
    chk('dim %d claim matches RESULTS.md' % d, tab4[d][1] == results[d][1],
        '%d vs %d' % (tab4[d][1], results[d][1]))
    fv, _ = floor_for(d)
    chk('dim %d prev matches published.py' % d, tab4[d][0] == fv,
        'tex %d vs floor_for %d' % (tab4[d][0], fv))
    chk('dim %d claim beats prev' % d, tab4[d][1] > tab4[d][0],
        '%d <= %d' % (tab4[d][1], tab4[d][0]))

# monotonicity of the claims taken together with the floors
best = {}
for d in range(1, 97):                     # 97, not 96: dimension 96 is now a claim
    try:
        best[d] = floor_for(d)[0]
    except ValueError:
        pass
for d, (_, v, _) in tab4.items():
    best[d] = max(best[d], v)
for d in range(2, 97):
    if d in best and d - 1 in best:
        chk('monotone at %d' % d, best[d] >= best[d - 1],
            'K(%d)=%d < K(%d)=%d' % (d, best[d], d - 1, best[d - 1]))

# ------------------------------------------------------------------------- Table 1
print('== Table 1 (summary) ==')
t1 = TEX[TEX.index('\\label{tab:summary}'):
         TEX.index('\\end{tabular}', TEX.index('\\label{tab:summary}'))]
cap = TEX[TEX.index('\\caption{The forty-seven improvements'):TEX.index('\\label{tab:summary}')]
GROUPS = [
    ('$\\Gsz$ cross-sections', [68, 69, 70, 71]),
    ('layers over $\\Leech$', [25, 27, 38]),
    ('layers over $\\Pff$', list(range(49, 62))),
    ('the \\textsc{ers} chain', [62, 63]),
    # Dimension 96 is the same construction at a different chain, (96,24,6,1) rather
    # than (n,15,2), and it gets its OWN row.  Merged into the row above, the factor
    # extremes would fall at 62 and 63 while the dimension range ran to 96, breaking the
    # caption's claim that both extremes sit at an end of the range -- the check below.
    ('the \\textsc{ers} chain at $96$', [96]),
    ('layers over $\\Gsz$', list(range(73, 96))),
    ('constant-weight codes', [39]),
]
allg = sorted(sum((d for _, d in GROUPS), []))
chk('Table 1 groups partition Table 4', allg == sorted(tab4),
    'groups %d dims, table4 %d' % (len(allg), len(tab4)))
chk('Table 1 total row says 47', 'total & & $47$' in t1)
smallest = []
for name, dims in GROUPS:
    row = [l for l in t1.split('\n') if l.startswith(name + ' &')]
    chk('Table 1 has a row for %r' % name, len(row) == 1, '%d rows' % len(row))
    if len(row) != 1:
        continue
    cells = row[0].split('&')
    chk('%s count column' % name, num(cells[2]) == len(dims),
        'says %s, is %d' % (cells[2], len(dims)))
    gains = sorted(results[d][1] - results[d][0] for d in dims)
    facs = sorted(results[d][1] / results[d][0] for d in dims)
    smallest.append(gains[0])
    chk('%s largest gain' % name, num(cells[3]) == gains[-1],
        '%s vs %d' % (cells[3].strip(), gains[-1]))
    fm = re.findall(r'(\d+\.\d+)', cells[4])
    chk('%s factor lo' % name, abs(float(fm[0]) - facs[0]) < 0.005,
        '%s vs %.4f' % (fm[0], facs[0]))
    chk('%s factor hi' % name, abs(float(fm[-1]) - facs[-1]) < 0.005,
        '%s vs %.4f' % (fm[-1], facs[-1]))
    if len(dims) > 1:
        # extremes must sit at the ends of the dimension range, as the caption asserts
        ends = {dims[0], dims[-1]}
        gmax = max(dims, key=lambda d: results[d][1] - results[d][0])
        fmin = min(dims, key=lambda d: results[d][1] / results[d][0])
        fmax = max(dims, key=lambda d: results[d][1] / results[d][0])
        # the caption's claim is about the DISPLAYED two-decimal factors: a dimension whose
        # factor rounds to the same value as an end of the range is not a counterexample
        # (dimensions 25 and 27 both show 1.00 since 2026-09-07)
        _r = lambda d: round(results[d][1] / results[d][0], 2)
        chk('%s extremes at the ends of the range' % name,
            gmax in ends and _r(fmin) == min(_r(dims[0]), _r(dims[-1]))
            and _r(fmax) == max(_r(dims[0]), _r(dims[-1])),
            'gmax %d fmin %d fmax %d' % (gmax, fmin, fmax))
chk('the caption carries no figures of its own',
    not re.findall(r'\$\+[\d\\,]+\$', cap), 'a numeral crept back into the caption')
# The prose under Table 1 counts the rows whose factor exceeds 2.  It said "Two" until
# dimension 96 arrived at 2.07 and made it three -- a counted claim about a table, sitting
# directly under that table, which nothing checked.  Count them from the table itself.
_over2 = [name for name, dims in GROUPS
          if max(results[d][1] / results[d][0] for d in dims) > 2]
_WORD = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five', 6: 'Six'}
# The companion sentence claimed the near-1 rows also carried the largest absolute gains.
# They do not, and did not before dimension 96 either: they top out at 254 553 142 while the
# cross-sections row reaches 2 271 920 766.  Check both halves of the replacement.
_near1 = [(name, dims) for name, dims in GROUPS
          if max(results[d][1] / results[d][0] for d in dims) < 1.2]
_span = sum(len(dims) for _n, dims in _near1)
chk('the near-1 rows span the stated number of dimensions',
    '$%d$ of the $%d$' % (_span, len(tab4)) in TEX,
    'they span %d of %d' % (_span, len(tab4)))
_bygain = sorted(GROUPS, reverse=True,
                 key=lambda nd: max(results[d][1] - results[d][0] for d in nd[1]))
_top2over2 = all(max(results[d][1] / results[d][0] for d in dims) > 2
                 for _n, dims in _bygain[:2])
chk('the two largest absolute gains are in rows with factor above 2', _top2over2,
    'top two by gain: %s' % ', '.join(n for n, _d in _bygain[:2]))

chk('the prose counts the rows with factor above 2',
    '%s rows exceed $2$' % _WORD[len(_over2)] in TEX,
    'there are %d such rows (%s), and the paper does not say so'
    % (len(_over2), ', '.join(_over2)))

chk('the caption states where the extremes fall',
    'both extremes attained at an end of the dimension range' in FLAT)

# ------------------------------------------------------- one-point distributions
print('== Proposition 4.8 (one-point distributions) ==')
ONE = {
    'E_8':    (8, 2, 240, [126, 56]),
    'Leech':  (24, 4, 196560, [93150, 47104, 4600]),
    'P48':    (48, 6, 52416000, [23766960, 12608784, 1678887, 36848]),
    'Gamma72': (72, 8, 6218175600, [2603658750, 1512243200, 280928256, 13959168, 127800]),
}
for L, (n, mu, N, row) in ONE.items():
    chk('%s one-point has mu/2+1 entries' % L, len(row) == mu // 2 + 1,
        '%d vs %d' % (len(row), mu // 2 + 1))
    tot = row[0] + 2 * sum(row[1:]) + 2
    chk('%s one-point sums to N' % L, tot == N, '%d vs %d' % (tot, N))
    # every moment identity, degrees 0,2,4,6,8,10 (7-design for E8 so only to 6)
    top = 11 if n % 24 == 0 else 7
    for p in range(0, top + 1, 2):
        lhs = sum((F(c, mu) ** p) * (row[c] * (2 if c else 1)) for c in range(len(row)))
        lhs += 2 * F(1, 1) ** p  # the two vectors +-u themselves, cos = +-1
        rhs = F(N) * F(1)
        den = 1
        for j in range(p // 2):
            den *= (n + 2 * j)
        # E[<u,Z>^p] = (p-1)!! |u|^p ; normalised u so = (p-1)!!
        dfac = 1
        for j in range(1, p, 2):
            dfac *= j
        rhs = F(N * dfac, den)
        chk('%s moment p=%d' % (L, p), lhs == rhs, '%s vs %s' % (lhs, rhs))

# ------------------------------------------------------------------ cap arithmetic
print('== Cap construction arithmetic ==')
chk('dim 25 level-3/4 count', 196560 + 2 * 248 * (2 - 1) + 2 == 197058)
chk('dim 25 claim is the lens-head configuration', 196560 - 1006 + 1 + 2 * 1006 + 2 == 197569)
chk('paper says 197,569', '197\\,569' in TEX)
chk('dim 25 with a 249-line class', 196560 + 2 * 249 * (2 - 1) + 2 == 197060)
chk('paper says 197060', '197\\,060' in TEX)
chk('dim 27 count', 196560 + 2 * 4 * 248 * (3 - 1) + 12 == 200540)
chk('dim 27 arithmetic pieces', 2 * 4 * 248 * 2 == 3968 and 196560 + 3968 + 12 == 200540)
chk('dim 27 five classes delete 2480', 5 * 496 == 2480)
chk('dim 27 four classes delete 1984', 4 * 496 == 1984)
chk('dim 27 difference is 496', 2480 - 1984 == 496)
chk('dim 27 ceiling floor(2*12/3)=8', 2 * 12 // 3 == 8)
chk('dim 27 gain from 8 classes-worth', 196560 + 2 * 248 * 8 + 12 == 200540)
chk('dim 38 equator empty', 196560 - 2 * 98280 + 2 * 98280 * 3 == 589680)
chk('dim 38 ceiling is 3K(24)', 6 * 98280 == 3 * 196560 == 589680)
chk('dim 38 total', 589680 + 1932 == 591612)
chk('dim 38 arena size', 6120 < 98280)
chk('Leech lines', 196560 // 2 == 98280)
chk('C(200540,2)', comb(200540, 2) == 20108045530, str(comb(200540, 2)))
chk('C(3324,2)', comb(3324, 2) == 5522826, str(comb(3324, 2)))

# floor(K(k)/3) for k=1..7
kk = [COHN[k] for k in range(1, 8)]
chk('K(1..7)', kk == [2, 6, 12, 24, 40, 72, 126], str(kk))
chk('floor(K(k)/3) list', [k // 3 for k in kk] == [0, 2, 4, 8, 13, 24, 42],
    str([k // 3 for k in kk]))

# calibration table
print('== Table 3 (calibration) ==')
CAL = {1: ('1p', 197058), 2: ('2t', 198550), 3: ('4t', 200540), 4: ('8t', 204520),   # k=1: the level-3/4 construction, not the claim
       5: ('12t,2p', 209496), 6: ('24t', 220440), 7: ('42t', 238350)}
for k, (parts, val) in CAL.items():
    nt = int(re.search(r'(\d+)t', parts).group(1)) if 't' in parts else 0
    npr = int(re.search(r'(\d+)p', parts).group(1)) if 'p' in parts else 0
    W = 3 * nt + 2 * npr
    chk('calib k=%d |W| == K(%d)' % (k, k), W == COHN[k], '%d vs %d' % (W, COHN[k]))
    got = 196560 + 2 * 248 * (nt * 2 + npr * 1) + W
    chk('calib k=%d count' % k, got == val, '%d vs %d' % (got, val))
    if k != 1:
        chk('calib k=%d vs published' % k,
            (got == COHN[24 + k]) == (k not in (1, 3)),
            'count %d, published %d' % (got, COHN[24 + k]))
chk('calib diffs +2 and +496', 197058 - COHN[25] == 2 and 200540 - COHN[27] == 496,
    '%d %d' % (197058 - COHN[25], 200540 - COHN[27]))

# ------------------------------------------------------------ Corollary 4.3 algebra
print('== Corollary 5.2 (pairs versus triples) ==')
for gamma, tmin, cmax, mmax in [(F(1, 4), F(2, 3), F(-1, 2), 3), (F(1, 3), F(3, 4), F(-1), 2)]:
    chk('t_min at gamma=%s' % gamma, F(1, 2) / (1 - gamma) == tmin,
        '%s vs %s' % (F(1, 2) / (1 - gamma), tmin))
    got = (F(1, 2) - tmin) / (1 - tmin)
    chk('cos bound at gamma=%s' % gamma, got == cmax, '%s vs %s' % (got, cmax))
    chk('paper formula gamma/(2gamma-1) at gamma=%s' % gamma,
        gamma / (2 * gamma - 1) == cmax and 'gamma/(2\\gamma-1)'.replace('gamma','\\gamma')
        or True, '')
    chk('tex carries gamma/(2gamma-1) at gamma=%s' % gamma,
        TEX.count(r'\gamma/(2\gamma-1)') == 2, 'count %d' % TEX.count(r'\gamma/(2\gamma-1)'))
    chk('correct formula gamma/(2gamma-1) at gamma=%s' % gamma,
        gamma / (2 * gamma - 1) == cmax, '%s' % (gamma / (2 * gamma - 1)))
    chk('m <= 1-1/c at gamma=%s' % gamma, floor(1 - 1 / cmax) == mmax,
        '%s' % (1 - 1 / cmax))
chk('gamma = floor(mu/3)/mu for Leech', F(4 // 3, 4) == F(1, 4))
chk('gamma = floor(mu/3)/mu for P48', F(6 // 3, 6) == F(1, 3), '%s' % F(6 // 3, 6))
chk('the old formula floor(mu/4)/mu is gone', 'mu/4' not in TEX.replace('\\', ''),
    'floor(mu/4) still in the tex')
chk('gamma = floor(mu/3)/mu for Gamma72', F(8 // 3, 8) == F(1, 4))
chk('gamma=3/8 would kill every pair', (F(1, 2) - F(4, 5)) / (1 - F(4, 5)) == F(-3, 2),
    '%s' % ((F(1, 2) - F(4, 5)) / (1 - F(4, 5))))
chk('gamma=3/8 forces t>=4/5', F(1, 2) / (1 - F(3, 8)) == F(4, 5))
chk('b at t=3/4 is 1/2', F(1, 4) == (1 - F(3, 4)))
chk('c6 automatic at t=3/4', sqrt(1 - 0.75) * 1.0 <= 0.5 + 1e-12)
chk('c6 fails at t=2/3 for coincident directions', sqrt(1 - 2 / 3) * 1.0 > 0.5)
chk('30 degrees at t=2/3', abs(0.5 / sqrt(1 - 2 / 3) - sqrt(3) / 2) < 1e-12)
chk('t=3/4 keeps Leech class: 3/4*1/4+1/4 <= 1/2',
    F(3, 4) * F(1, 4) + F(1, 4) <= F(1, 2), '%s' % (F(3, 4) * F(1, 4) + F(1, 4)))
chk('7/16 is that value', F(3, 4) * F(1, 4) + F(1, 4) == F(7, 16))

# ----------------------------------------------------------------- Edel-Rains-Sloane
print('== Edel-Rains-Sloane arithmetic ==')
chk('ERS dim 64 reverse-engineers',
    2 ** 28 + 30828 * 2048 + 10416 * 16 + 128 == 331737984,
    str(2 ** 28 + 30828 * 2048 + 10416 * 16 + 128))
chk('ERS dim 80 reverse-engineers',
    2 ** 30 + 143780 * 2048 + 20540 * 16 + 160 == 1368532064,
    str(2 ** 30 + 143780 * 2048 + 20540 * 16 + 160))
chk('dim 39 chain', 327680 + 128 * 3324 + 4 * comb(39, 2) == 756116,
    str(327680 + 128 * 3324 + 4 * comb(39, 2)))
chk('dim 39 pieces', 128 * 3324 == 425472 and 4 * 741 == 2964 and comb(39, 2) == 741)
chk('dim 39 old record', 327680 + 128 * 3323 + 2964 == COHN[39],
    '%d vs %d' % (327680 + 128 * 3323 + 2964, COHN[39]))
chk('dim 39 gain is 128', 756116 - COHN[39] == 128, str(756116 - COHN[39]))
chk('ceil(39/4)=10', -(-39 // 4) == 10)
chk('A(8,2)=128 is 2^7', 2 ** 7 == 128)
chk('A(2,1)=4', 2 ** 2 == 4)
chk('2^32 vs Gamma72', 2 ** 32 < 6218175600,
    '2^32 = %d is BELOW 6218175600 -- check the stopping argument' % 2 ** 32)
chk('2^30 far below Gamma72', 2 ** 30 < 6218175600)

# ------------------------------------------------------------------- classes / Caro-Wei
print('== Proposition 5.6 (classes) ==')
LINES = {'Leech': 98280, 'P48': 26208000, 'Gamma72': 3109087800}
DEG = {'Leech': 4600, 'P48': 36848, 'Gamma72': 13959168 + 127800}
CW = {'Leech': 22, 'P48': 712, 'Gamma72': 221}
_G72_DATA = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
_P48_DATA = os.path.join(KV, 'verifications', 'improved', 'dim49-63-p48-caps', 'data')


def _largest_class(*paths):
    """The largest independent set the package ships, over the size tables named.

    Written down here as 2106 until 2026-08-27, which was the best the orbit-pool search
    found; the greedy family rebuilt the day before holds one of 2118, so the paper had gone
    stale against its own repository and this guard could not see it.  A literal here agrees
    with the literal there whatever the shipped data says, which is the whole failure mode.
    """
    best = 0
    for p in paths:
        if os.path.exists(p):
            best = max(best, max(_npy_ints(p)))
    return best


FOUND = {'Leech': 248,
         'P48': _largest_class(os.path.join(_P48_DATA, 'class_sizes_p48.npy')),
         'Gamma72': _largest_class(os.path.join(_G72_DATA, 'class_sizes_32000.npy'),
                                   os.path.join(_G72_DATA, 'class_sizes_gpu.npy'))}
RATIO = {'Leech': 11.3, 'P48': 9.93, 'Gamma72': 9.58}
chk('the shipped class-size tables are readable',
    FOUND['P48'] > 0 and FOUND['Gamma72'] > 0,
    'no size table under %s or %s' % (_P48_DATA, _G72_DATA))
chk('the paper quotes the largest independent set each family holds',
    all(_grp(FOUND[_L]) in TEX for _L in ('P48', 'Gamma72')),
    'measured %d and %d' % (FOUND['P48'], FOUND['Gamma72']))
for L in LINES:
    chk('%s line count' % L, LINES[L] * 2 == ONE[L if L != 'Leech' else 'Leech'][2],
        '%d' % LINES[L])
    cw = -(-LINES[L] // (DEG[L] + 1))
    chk('%s Caro-Wei' % L, cw == CW[L], 'computed %d, paper %d' % (cw, CW[L]))
    r = FOUND[L] / CW[L]
    chk('%s ratio' % L, abs(round(r,1) - RATIO[L]) < 1e-9 or abs(round(r,2)-RATIO[L])<1e-9, 'computed %.3f, paper %s' % (r, RATIO[L]))
chk('P48 degree is n[3]', DEG['P48'] == 36848)
_p48 = os.path.join(KV, 'verifications', 'improved', 'dim49-63-p48-caps',
                    'data', 'class_sizes_p48.npy')
if os.path.exists(_p48):
    _s = _npy_ints(_p48)
    _r, _l = len(_s), sum(_s)
    chk('P_48 family size matches class_sizes_p48.npy',
        _grp(_r) in TEX and _grp(_l) in TEX, '%d sets, %d lines' % (_r, _l))
    chk('P_48 coverage percentage', '%.1f' % (100.0 * _l / 26208000) in TEX,
        '%.1f' % (100.0 * _l / 26208000))
else:
    chk('class_sizes_p48.npy present', False, 'missing: %s' % _p48)

# ---------------------------------------------------------------- the two Gamma_72 families
# Section 5.6 describes both disjoint families over Gamma_72 by size, coverage and the range
# of set sizes.  Those five numbers went stale when the greedy family was extended from
# 31 000 sets to 32 000 on 2026-08-26: the paper still said 31 000 sets, 63 464 272 lines and
# sizes 1 938 to 2 111 against the shipped 32 000 / 65 570 602 / 1 932 to 2 118.  No BOUND
# moved -- final81-95.py reads the .npy files and never the prose -- but the paper described
# an object the repository no longer contained.  Anchored on the sentences rather than on
# "appears somewhere in the file": both families now have 32 000 sets, so a bare membership
# test would pass while reading the wrong one.
_G72 = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
_img_p = os.path.join(_G72, 'class_sizes_32000.npy')
_gpu_p = os.path.join(_G72, 'class_sizes_gpu.npy')
if os.path.exists(_img_p) and os.path.exists(_gpu_p):
    _img, _gpu = sorted(_npy_ints(_img_p)), sorted(_npy_ints(_gpu_p))
    _flat = re.sub(r'[~%ss]+' % BS, ' ', TEX)                # LaTeX ties and line breaks
    _E = re.escape                                           # BS+'Gsz' is not a regex escape
    _NUM = r'[$]([^$]+)[$]'

    def _sent(pat, what):
        m = re.search(pat, _flat)
        chk('section 5.6 states the %s' % what, m is not None, 'sentence not found')
        return m

    _m = _sent(_E('and $%s$ of $%sGsz$, covering ' % (_grp(len(_img)), BS)) + _NUM
               + _E(' lines'), 'image family of Gamma_72')
    if _m:
        chk('image family coverage matches class_sizes_32000.npy',
            _m.group(1) == _grp(sum(_img)),
            'paper %s, file %s' % (_m.group(1), _grp(sum(_img))))
    _m = _sent(_E('over $%sGsz$ it gives ' % BS) + _NUM + _E(' sets covering ') + _NUM
               + _E(' lines'), 'greedy family of Gamma_72')
    if _m:
        chk('greedy family size matches class_sizes_gpu.npy',
            (_m.group(1), _m.group(2)) == (_grp(len(_gpu)), _grp(sum(_gpu))),
            'paper %s sets / %s lines, file %s / %s'
            % (_m.group(1), _m.group(2), _grp(len(_gpu)), _grp(sum(_gpu))))
    _m = _sent(_E('begin at ') + _NUM + _E(' lines and decay to ') + _NUM
               + _E(', while the sets found afresh lie between ') + _NUM + _E(' and ')
               + _NUM, 'size range of both families')
    if _m:
        chk('image family size range matches the file',
            (_m.group(1), _m.group(2)) == (_grp(_img[-1]), _grp(_img[0])),
            'paper %s..%s, file %s..%s' % (_m.group(2), _m.group(1),
                                           _grp(_img[0]), _grp(_img[-1])))
        chk('greedy family size range matches the file',
            (_m.group(3), _m.group(4)) == (_grp(_gpu[0]), _grp(_gpu[-1])),
            'paper %s..%s, file %s..%s' % (_m.group(3), _m.group(4),
                                           _grp(_gpu[0]), _grp(_gpu[-1])))
    # and the crossover the paper names: the first dimension the greedy family wins
    _rows = os.path.join(_G72, '..', 'rows81-95.json')
    if os.path.exists(_rows):
        _rw = json.load(io.open(_rows, encoding='utf-8'))
        _ci, _cg = [0], [0]
        for _v in sorted(_img, reverse=True):
            _ci.append(_ci[-1] + _v)
        for _v in sorted(_gpu, reverse=True):
            _cg.append(_cg[-1] + _v)
        _over = [r['dim'] for r in _rw
                 if r['T'] + r['P'] <= min(len(_img), len(_gpu))
                 and _cg[r['T'] + r['P']] > _ci[r['T'] + r['P']]]
        _m = _sent(_E('the crossover falling at dimension ') + _NUM, 'crossover dimension')
        if _m and _over:
            chk('the crossover dimension is where the greedy family overtakes',
                _m.group(1) == _grp(min(_over)),
                'paper %s, measured %d' % (_m.group(1), min(_over)))
else:
    chk('the Gamma_72 class-size files are present', False, 'missing under %s' % _G72)

# ---------------------------------------------------------------- Table 2 vs validate_lp
# Table 2 is the paper's validation of the cross-section programme.  It carries ONE MORE ROW
# than validate_lp.py has cases, and that is deliberate: the Leech k = 3 orthogonal Gram
# mu*I_3 admits two different embeddings with two different counts.  In coordinates where a
# minimal vector has raw norm 32, 4e1+4e2, 4e3+4e4, 4e5+4e6 gives 19530 while
# 4e1+4e2, 4e1-4e2, 4e3+4e4 gives 19962; both triples are pairwise orthogonal of norm mu, so
# both have Gram 32*I_3.  Verified 2026-08-26 by direct count over all 196560 minimal vectors
# from an independent Leech construction.  **The cross-section count is therefore NOT a
# function of the Gram matrix** -- which is exactly why the k-point LP, which sees only the
# Gram, brackets every embedding rather than pinning one, and needs a realizability
# argument to be attained (section 4.5).  validate_lp.py records the 19530 embedding; the paper records both.
# An earlier version of this file asserted that 19962 "was never a count" and that the row
# should be removed.  That was wrong, and validate_lp.py's own note says so.
print('== Table 2 (validation) against common/validate_lp.py ==')
_B = chr(92)
_vl = io.open(os.path.join(KV, 'common', 'validate_lp.py'), encoding='utf-8').read()
_cases = _vl[_vl.index('CASES = ['):_vl.index('nbelow = nabove = 0')]
_truths = sorted(int(m) for m in re.findall(
    ',' + _B + 's*(' + _B + 'd+),' + _B + 's*(?:7|11),', _cases))
_i = TEX.index('label{tab:validate}')
_t2 = TEX[_i:TEX.index('end{tabular}', _i)]
_paper = sorted(num(m) for m in re.findall('mathbf{([^}]+)}', _t2))
_extra = [v for v in _paper if v not in _truths]
chk('Table 2 covers every validate_lp case', [v for v in _truths if v not in _paper] == [],
    'validate_lp cases missing from the paper: %s'
    % [v for v in _truths if v not in _paper])
chk('the only extra Table 2 row is the second orthogonal embedding', _extra == [19962],
    'paper carries %s beyond validate_lp' % _extra)
chk('both orthogonal counts are present', 19530 in _paper and 19962 in _paper,
    'the mu*I_3 Gram has two embeddings, 19530 and 19962; Table 2 records both')

# ------------------------------------------------------------------------ obstructions
print('== Obstructions ==')
chk('dim 19: 10668+4*320', 10668 + 4 * 320 == 11948 == COHN[19])
chk('dim 20/21 128*16', 128 * 16 == 2048)
chk('dim 17: 5346+2*192', 5346 + 2 * 192 == 5730 == COHN[17],
    '%d vs COHN %d' % (5346 + 2 * 192, COHN[17]))
chk('dim 17: alpha = 192 is exact, so 5730 is the family maximum',
    5346 + 2 * 192 == 5730 and 13 + 5 + 45 + 56 + 45 + 27 == 191 and 1 + 191 == 192)
chk('dim 17 +384', 2 * 192 == 384)
chk('class upper bound 4680/11', floor(F(4680, 11)) == 425, str(F(4680, 11)))
chk('4680/11 decimal', abs(4680 / 11 - 425.4545) < 1e-3)
chk('Hall 1024/2025 + 44/891 < 1', F(1024, 2025) + F(44, 891) < 1,
    '%s = %.4f' % (F(1024, 2025) + F(44, 891), float(F(1024, 2025) + F(44, 891))))
chk('Lambda_23 kissing is 93150', COHN[23] == 93150)
# raising the class record by one line: gain per dimension
GAINS = {25: 2, 26: 8, 27: 16, 28: 32, 29: 52, 30: 96, 31: 168}
chk('paper lists the corrected gains',
    '$+2$, $+8$, $+16$, $+32$, $+52$, $+96$ and $+168$' in FLAT)
for k in range(1, 8):
    d = 24 + k
    # gain = 2 * (extra line) * sum_i (|Z_i| - 1) over the classes used
    parts = CAL[k][0]
    nt = int(re.search(r'(\d+)t', parts).group(1)) if 't' in parts else 0
    npr = int(re.search(r'(\d+)p', parts).group(1)) if 'p' in parts else 0
    g = 2 * (nt * 2 + npr * 1)
    chk('one extra class line gives +%d in dim %d' % (GAINS[d], d), g == GAINS[d],
        'computed %d, paper %d' % (g, GAINS[d]))

# ------------------------------------------------------------------------- Ozeki
print('== Ozeki cross-checks ==')
chk('Ozeki degree-3 ratio is 12309600',
    F(23775066324172800000, 1931424768000) == 12309600,
    '%s' % F(23775066324172800000, 1931424768000))
chk('458386320/43 is not an integer', F(458386320, 43).denominator != 1,
    '%s' % F(458386320, 43))
chk('Ozeki line: low end exact', 6732912 - 164 * 278 == 6687320)
chk('Ozeki line: high end is 6702736', 6732912 - 164 * 184 == 6702736)
chk('paper states 6702736', '6\\,702\\,736' in TEX)
chk('the stale certified 6702737 is gone', '6\\,702\\,737' not in TEX,
    'common/validate_lp.py now certifies 6702736 exactly, on the vertex')

# --------------------------------------------------- dimension 38 sits AT its ceiling
# The pole layer read 1908 until 2026-08-27, when a rotation clearing all 1932 was found by
# descent rather than by sampling.  The paper now asserts something stronger than a number:
# that the total IS the family ceiling 3K(24)+K(14).  Guard both halves, and guard that the
# superseded pair 589680+1908 has not survived anywhere in the tex.
print('== Dimension 38 at the ceiling ==')
chk('dim 38 ceiling identity', 589680 + 1932 == 591612 and 589680 == 3 * 196560
    and 1932 == 1932)
chk('the paper says all 1932 rotated points are poles',
    re.search(r'all@n?@s*[$]1@,932[$] of them are poles'.replace('@n', chr(92) + 'n')
              .replace('@s', chr(92) + 's').replace('@,', chr(92) + chr(92) + ','),
              FLAT) is not None,
    'the sentence that carries the +24 is gone')
chk('the paper states the ceiling 3K(24)+K(14)',
    '3K(24)+K(14)' in FLAT.replace(' ', ''),
    'the claim is that the construction is AT its ceiling, not near it')
chk('the superseded 589680+1908 is gone from the tex',
    '589' + chr(92) + ',680+1' + chr(92) + ',908' not in TEX)
chk('the superseded 591588 is gone from the tex',
    '591' + chr(92) + ',588' not in TEX)
# The three numbers in the sampling sentence have to agree with ONE cap measure, and until
# 2026-08-27 they did not: the paper read "1.24% of S^13 ... loses 56 points on average ...
# gives e^{-24}", and 1.24% is 24/1932 -- back-derived from the single rotation that lost 24,
# not computed.  A 30-degree cap on S^13 has relative measure
# Int_0^{pi/6} sin^12 / Int_0^pi sin^12, which is 1.4978e-5, so the 1932 caps cover 2.89% and
# the expected loss is 55.9, matching the MEASURED 56.  Derive the measure here from
# sin^12 alone and require all three of the paper's figures to follow from it.
_S13 = 0.0
_N = 240000
for _i in range(_N + 1):                      # Simpson on sin^12 over [0, pi/6] and [0, pi]
    _x = math.pi * _i / _N
    _w = 1 if _i in (0, _N) else (4 if _i % 2 else 2)
    _S13 += _w * math.sin(_x) ** 12
_num = 0.0
for _i in range(_N + 1):
    _x = (math.pi / 6) * _i / _N
    _w = 1 if _i in (0, _N) else (4 if _i % 2 else 2)
    _num += _w * math.sin(_x) ** 12
_capfrac = (_num / 6) / _S13                  # the step lengths cancel except for the 1/6
_m = re.search(r'relative measure\s+[$]([\d.]+)' + BS + BS + r'times10\^[{]-5[}][$]', FLAT)
chk('the paper states the 30-degree cap measure on S^13',
    bool(_m) and abs(float(_m.group(1)) * 1e-5 - _capfrac) < 5e-8,
    'paper %s, computed %.4g' % (_m.group(1) if _m else None, _capfrac * 1e5))
_m = re.search(r'[$]1' + BS + BS + r',932[$] caps cover at most [$]([\d.]+)' + BS + BS
               + r'%[$]', FLAT)
chk('and the coverage that follows from it',
    bool(_m) and abs(float(_m.group(1)) / 100 - 1932 * _capfrac) < 5e-5,
    'paper %s%%, computed %.4f%%' % (_m.group(1) if _m else None, 100 * 1932 * _capfrac))
_m = re.search(r'loses [$](\d+)[$] points on average', FLAT)
chk('and the measured average loss agrees with that coverage',
    bool(_m) and abs(int(_m.group(1)) - 1932 * 1932 * _capfrac) < 1.0,
    'paper %s, computed %.1f' % (_m.group(1) if _m else None, 1932 * 1932 * _capfrac))
_m = re.search(r'e\^[{]-(\d+)[}]', TEX)
chk('and the independent-points exponent is that loss, not the old 24',
    bool(_m) and abs(int(_m.group(1)) - 1932 * 1932 * _capfrac) < 1.0,
    'paper e^{-%s}, computed e^{-%.1f}'
    % (_m.group(1) if _m else None, 1932 * 1932 * _capfrac))
chk('the paper exhibits the Cayley transform it used, not just its name',
    '(I-S)(I+S)^{-1}' in FLAT.replace(' ', ''))
chk('the paper states the exact integer test the certificate passes',
    'isqrt(48D^2)' in FLAT.replace(' ', '').replace(chr(92) + 'operatorname', '')
    .replace('{', '').replace('}', '').replace(chr(92) + ',', ''),
    'the threshold is an integer square root, which is what makes it exact')

# ------------------------------------------- the axis layers of 75-95 are rotated copies
# They were lattice-shell sets until 2026-08-27 -- 248 points at k = 23 against tau(23) =
# 93 150 -- and are now rotated copies of the cap directions, worth +213 230 spheres over
# nineteen dimensions.  Guard the sentences that carry that, and that the superseded claim
# that a rotation does not produce a pole set is gone.
print('== The rotated axis layers ==')
chk('the paper says a rotated copy is a 60-degree code for free',
    'rotated copy of the direction set is automatically' in FLAT)
# The sentence that carries this is one calculation, not three literals, and the version
# shipped before 2026-08-27 got every noun wrong: "a 30-degree cap covers 2.3e-3 of S^22,
# so a rotation of the 93 150 directions loses about two hundred".  2.3e-3 is the UNION of
# the 93 150 caps, five orders of magnitude off one cap, and "two hundred" is the
# prediction rather than the loss -- 93 150 - 200 is not 93 044.  So derive all three.
#
#   fraction of S^{n-1} within 30 degrees of a point = int_0^{pi/6} sin^{n-2} / int_0^pi.
# For n = 23 the exponent 21 is odd, so u = cos t turns both into polynomials and the
# whole thing is exact in Q + Q*sqrt3: with p = 10,
#   int_c^1 (1-u^2)^p du = A - B*sqrt3/2,  A = sum_j C(p,j)(-1)^j/(2j+1),
#                                          B = sum_j C(p,j)(-1)^j (3/4)^j/(2j+1),
# and int_{-1}^1 = 2A.
_p = (23 - 2 - 1) // 2
_A = sum(F(comb(_p, j) * (-1) ** j, 2 * j + 1) for j in range(_p + 1))
_B = sum(F(comb(_p, j) * (-1) ** j, 2 * j + 1) * F(3, 4) ** j for j in range(_p + 1))
_cap = float(_A - _B * F(sqrt(3)).limit_denominator(10 ** 12) / 2) / float(2 * _A)
_nz = 93150
_union = _nz * _cap
_eloss = _nz * _union
_layer = _shipped_poles(23)          # from the package's own data, not from the paper
chk('a single 30-degree cap on S^22 is 2.3e-8, and the paper says so',
    abs(_cap - 2.284e-8) < 1e-11 and '2.3' + BS + 'times10^{-8}' in TEX,
    'measured %.4g' % _cap)
chk('the 93 150 caps together cover at most 2.1e-3, and the paper says so',
    abs(_union - 2.128e-3) < 1e-5 and '2.1' + BS + 'times10^{-3}' in TEX,
    'measured %.4g' % _union)
chk('the expected loss of a Haar rotation is 198, and the paper says so',
    197.0 < _eloss < 199.0 and '$198$' in TEX,
    'measured %.1f' % _eloss)
chk('the paper does not call the union of the caps a single cap',
    '2.3' + BS + 'times10^{-3}' not in TEX,
    'the pre-2026-08-27 sentence, which named the wrong quantity')
# The sentence names three numbers and they are one calculation.  Read the loss out of the
# paper rather than writing it in here: a literal in a check makes today's value a
# requirement, which is exactly how the cap measure came to be un-correctable.
_m_loss = re.search(r'loses \$(\d+)\$, leaving \$([0-9' + re.escape(BS) + r',]+)\$', TEX)
chk('the paper states the loss and the layer it leaves',
    _m_loss is not None, 'the sentence has changed shape')
if _m_loss:
    _said_loss, _said_layer = int(_m_loss.group(1)), num(_m_loss.group(2))
    chk('the loss the paper states is the one the shipped layer implies',
        _said_loss == _nz - _layer,
        'the paper says %d, the data says %d' % (_said_loss, _nz - _layer))
    chk('the sentence closes: |Z| - loss = the layer',
        _nz - _said_loss == _said_layer and _said_layer == _layer,
        '%d - %d = %d against a shipped layer of %d'
        % (_nz, _said_loss, _nz - _said_loss, _layer))
    chk('the stated loss is below the bound it is derived from',
        _said_loss < _eloss, '%d lost against a bound of %.1f' % (_said_loss, _eloss))
# A bare % is a LaTeX comment: it silently deletes the rest of the line, and a script that
# rewrites a sentence containing $99\%$ can drop the backslash without anything looking
# wrong until the PDF is read.  That happened on 2026-08-29 and this is what found it.
_pct_bad = []
for _i, _ln in enumerate(TEX.replace(chr(13) + chr(10), chr(10)).split(chr(10))):
    if _ln.lstrip().startswith('%'):
        continue                                  # a whole-line comment is fine
    for _j, _ch in enumerate(_ln):
        if _ch == '%' and (_j == 0 or _ln[_j - 1] != BS):
            _pct_bad.append('line %d: %s' % (_i + 1, _ln.strip()[:60]))
            break
chk('no unescaped %% outside a comment line', not _pct_bad,
    '; '.join(_pct_bad[:3]))

chk('the paper names the Cayley transform it uses',
    '(I-S)(I+S)^{-1}' in FLAT.replace(' ', ''))
chk('the superseded sentence is gone',
    'an arbitrary rotation of the direction set' not in FLAT,
    'the paper used to say a rotation does not produce a pole set')
# Dimension 83 HAS a layer now -- all 604 of them, a rational rotation of a configuration
# that lives over Z[sqrt2] -- so the sentence this used to require is gone from the paper,
# and a guard for a sentence the paper no longer makes is a guard that can only mislead.
# What replaces it is the pair of checks above: the paper states 604 and the shipped layer
# is 604 = K(11).  Require the old claim to be ABSENT, so it cannot come back unnoticed.
chk('the paper no longer says a dimension is left without an axis layer',
    'Only $k=11$ is left' not in FLAT,
    'the superseded sentence about dimension 83 is still there')

# The size of the rotated layer, stated three ways in one sentence of S5.6.  It shipped as
# "between 99% and 100% of it" on 2026-08-27, which is false at seven of the eighteen k and
# thirteen points out at k = 8.  Recompute from the layers themselves.
_frac = _layer_fractions()
_pc = dict((k, 100.0 * p / z) for k, (z, p) in _frac.items())
chk('every dimension of 81-95 but 83 ships a rotated axis layer',
    len(_frac) >= 14, 'found %d' % len(_frac))
# The two percentages are READ from the paper.  Asserting them here is what makes a correct
# improvement fail the build -- the mistake this file made with the cap measure.
# The sentence used to carry a second, k>=12 threshold as well.  Every layer is now at
# 99.35% or better, so that clause said nothing the floor did not, and it is gone from the
# paper and from here together -- a guard for a sentence the paper no longer makes is a
# guard that passes on absence.
_mf = re.search(r'which is between \$(\d+)' + re.escape(BS) + r'%\$ and \$100'
                + re.escape(BS) + r'%\$ of it;', TEX)
chk('the paper states a floor for the size of the rotated layer', _mf is not None,
    'the fraction sentence has changed shape')
chk('the paper no longer carries the superseded k>=12 clause',
    'and above $' not in TEX.replace(chr(10), ' '),
    'the k>=12 threshold sentence is still there')
if _mf:
    _lo = int(_mf.group(1))
    _dlo = min(_pc.values())
    chk('the floor the paper states holds, and is not loose by more than a point',
        _lo <= _dlo < _lo + 1, 'paper %d%%, data %.2f%%' % (_lo, _dlo))
# ---- the record-configuration axis layers of dimensions 83, 91, 92 and 93 ---------------
# Four counts and three ceilings, all stated in the proof sketch of the Gamma_72 theorem.
# _shipped_poles() reads whichever layer file the dimension actually uses.
_rec = re.search(r'record configuration of \$.R\^k\$ instead gives\s+\$([0-9\\, ]+)\$, '
                 r'\$([0-9\\, ]+)\$ and \$([0-9\\, ]+)\$ poles', FLAT)
chk('the paper states the three record-source layers',
    _rec is not None, 'the sentence has changed shape')
if _rec:
    for _i, _k in enumerate((19, 20, 21)):
        _said = num(_rec.group(_i + 1))
        chk('dimension %d states its axis layer' % (72 + _k),
            _said == _shipped_poles(_k),
            'paper %d, shipped %d' % (_said, _shipped_poles(_k)))
_cei = re.search(r'against \$([0-9\\, ]+)\$, \$([0-9\\, ]+)\$ and \$([0-9\\, ]+)\$\.', FLAT)
chk('the paper states the three ceilings the layers are measured against',
    _cei is not None, 'the ceiling sentence has changed shape')
if _cei:
    for _i, _k in enumerate((19, 20, 21)):
        chk('the ceiling quoted for dimension %d is tau(%d)' % (72 + _k, _k),
            num(_cei.group(_i + 1)) == TAU[_k],
            'paper %d, table %d' % (num(_cei.group(_i + 1)), TAU[_k]))
_m11 = re.search(r'All \$([0-9\\, ]+) = K\(11\)\$ of them clear', FLAT)
chk('the paper states dimension 83 clears its whole configuration',
    _m11 is not None, 'the Z[sqrt2] sentence has changed shape')
if _m11:
    chk('dimension 83 states its axis layer, and it is K(11)',
        num(_m11.group(1)) == _shipped_poles(11) == TAU[11],
        'paper %d, shipped %d, tau %d' % (num(_m11.group(1)), _shipped_poles(11), TAU[11]))
# The paper says no configuration built in the standard R^k can be placed rationally in the
# axis spaces of k = 18 and 22.  That is a discriminant computation, and it is the reason
# those two dimensions are NOT in the list above -- so check the list has not grown.
chk('the paper names 18 and 22 as the axis spaces with no frame',
    'At $k=18$ and $k=22$ the discriminant is $3$' in FLAT,
    'the discriminant sentence has changed shape')
# ...and it must keep its SCOPE.  The sentence is about the standard metric; dimension 90's
# layer exists because the record configuration of R^18 is not written in that metric, so a
# reader who takes the sentence for a closure is reading it wrongly and the paper has to say
# which it means.  This is the "scoped claims lose their scope" failure, guarded.
chk('the paper says the no-frame sentence is about the metric, not the configuration',
    'That is a statement about the standard metric and not about the configuration' in FLAT,
    'the scope sentence after the discriminant claim is gone')

# ---- the axis layers of dimensions 85, 87, 89 and 90, from the PUBLISHED configurations --
# Four layers, four ceilings and four starting points in one sentence, plus the metrics that
# make the frames possible.  Every figure is read from the shipped files: the metric of each
# source is data/pole_src_<k>_c.npy, so the paper's diag(1^14,2) and the rest are checked
# against the array the search actually used and not against a transcription of it.
# k = 13 left this list on 2026-08-30: dimension 85's cap set became the Kappa
# section K_13, whose span is not Q-similar to the standard R^13 -- the Hasse
# invariants differ at 2 and 3 -- so no frame carries the record configuration into
# it and the layer is a rotated copy of K_13 instead.
_COHN_KS = (15, 17, 18)
_COHN_WAS = {15: 2340, 17: 5338, 18: 7374}


def _src_metric(k):
    # _npy_ints, not numpy: this file advertises itself as needing no third-party package
    p = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data',
                     'pole_src_%d_c.npy' % k)
    return [int(v) for v in _npy_ints(p)] if os.path.exists(p) else None


def _src_size(k):
    p = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data',
                     'pole_src_%d.npy' % k)
    return _npy_rows(p) if os.path.exists(p) else None


_num3 = r'[0-9' + re.escape(BS) + r', ]+'
_mlay = re.search(r'The layers\s+go from \$(' + _num3 + r')\$, \$(' + _num3 + r')\$ and \$('
                  + _num3 + r')\$ to \$(' + _num3 + r')\$, \$(' + _num3 + r')\$ and \$('
                  + _num3 + r')\$ poles, against \$K\(k\)\$ = \$(' + _num3 + r')\$, \$('
                  + _num3 + r')\$ and \$(' + _num3 + r')\$\.', FLAT)
chk('the paper states the three published-source layers, before and after, with their ceilings',
    _mlay is not None, 'the sentence has changed shape')
if _mlay:
    for _i, _k in enumerate(_COHN_KS):
        chk('dimension %d states the layer it started from' % (72 + _k),
            num(_mlay.group(_i + 1)) == _COHN_WAS[_k],
            'paper %d, was %d' % (num(_mlay.group(_i + 1)), _COHN_WAS[_k]))
        chk('dimension %d states the layer it ships' % (72 + _k),
            num(_mlay.group(_i + 4)) == _shipped_poles(_k),
            'paper %d, shipped %d' % (num(_mlay.group(_i + 4)), _shipped_poles(_k)))
        chk('the ceiling quoted for dimension %d is tau(%d)' % (72 + _k, _k),
            num(_mlay.group(_i + 7)) == TAU[_k],
            'paper %d, table %d' % (num(_mlay.group(_i + 7)), TAU[_k]))
        chk('dimension %d improved on what it started from' % (72 + _k),
            _shipped_poles(_k) > _COHN_WAS[_k],
            'shipped %d against %d' % (_shipped_poles(_k), _COHN_WAS[_k]))
# The three metrics, against the arrays the search used.  A metric is what decides whether a
# frame exists at all, so a wrong one in the paper would misstate the reason dimension 90
# works -- and nothing numeric elsewhere would catch it.
for _k, _pat in ((15, r'\operatorname{diag}(1^{14},2)'),
                 (17, r'\operatorname{diag}(1^{16},2)'),
                 (18, r'\operatorname{diag}(1^{16},2,6)')):
    _c = _src_metric(_k)
    _want = None
    if _c is not None:
        _tail = [v for v in _c if v != 1]
        _want = (r'\operatorname{diag}(1^{%d}%s)'
                 % (len(_c) - len(_tail), ''.join(',%d' % v for v in _tail)))
    chk('the metric the paper gives for k = %d is the one the source is written in' % _k,
        _c is not None and _want == _pat and _pat.replace(' ', '') in FLAT.replace(' ', ''),
        'shipped %s, paper %s' % (_want, _pat))
# ---- the Kappa cap sets, and the order-3 partitions -------------------------------------
# Dimension 85 stopped using the record configuration of R^13 on 2026-08-30, so the two
# guards that stood here -- how much of it is rational, and that the split closes to K(13) --
# ask about a claim the paper no longer makes.  These replace them, and they guard the claim
# that took its place: that the cap sets at k = 12 and 13 are the Kappa sections, that their
# triple counts are the shipped ones, and that k = 16 and 18 are now perfect.
def _capsize(k):
    p = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data',
                     'kap%d_W.npy' % k)
    return _npy_rows(p) if os.path.exists(p) else None


def _tri(fn):
    p = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data', fn)
    return _npz_rows(p, 'triples') if os.path.exists(p) else None


_mkap = re.search(r'The \$(' + _num3 + r')\$ minimal vectors of \$\\Leech\$ lying in the '
                  r'complement have\s+determinant \$3\^6\$ at minimal norm \$4\$', FLAT)
chk('the paper states the size of the K_12 section',
    _mkap is not None and num(_mkap.group(1)) == _capsize(12),
    'paper %s, shipped %s' % (num(_mkap.group(1)) if _mkap else None, _capsize(12)))
_mt12 = re.search(r'their partition into \$(' + _num3 + r')\$ zero-sum triples is the orbit '
                  r'set of \$\\pi\$', FLAT)
chk('the K_12 triple count in the paper is the shipped partition',
    _mt12 is not None and num(_mt12.group(1)) == _tri('kap12_parts.npz'),
    'paper %s, shipped %s' % (num(_mt12.group(1)) if _mt12 else None,
                              _tri('kap12_parts.npz')))
_mt13 = re.search(r'gives a \$13\$-dimensional section of \$(' + _num3 + r')\$', FLAT)
chk('the paper states the size of the K_13 section',
    _mt13 is not None and num(_mt13.group(1)) == _capsize(13),
    'paper %s, shipped %s' % (num(_mt13.group(1)) if _mt13 else None, _capsize(13)))
_mt13b = re.search(r'and those \$(' + _num3 + r')\$ triples are packed rather than given', FLAT)
chk('the K_13 triple count in the paper is the shipped partition',
    _mt13b is not None and num(_mt13b.group(1)) == _tri('kap13_parts.npz'),
    'paper %s, shipped %s' % (num(_mt13b.group(1)) if _mt13b else None,
                              _tri('kap13_parts.npz')))
_mom = re.search(r'from \$1\\,434\$ triples and \$7\$ pairs to \$(' + _num3 +
                 r')\$ triples, from \$2\\,457\$ and \$11\$ to \$(' + _num3 +
                 r')\$, from \$(' + _num3 + r')\$ and \$23\$ to \$(' + _num3 + r')\$, and from \$(' + _num3 +
                 r')\$ and \$49\$ to \$(' + _num3 + r')\$', FLAT)
chk('the paper states the four order-3 partitions, and they are the shipped ones',
    _mom is not None and num(_mom.group(1)) == _tri('parts_16.npz')
    and num(_mom.group(2)) == _tri('parts_18.npz')
    and num(_mom.group(4)) == _tri('parts_20.npz')
    and num(_mom.group(6)) == _tri('parts_22.npz'),
    'paper %s, shipped %s'
    % ((num(_mom.group(1)), num(_mom.group(2)), num(_mom.group(4)),
        num(_mom.group(6))) if _mom else None,
       (_tri('parts_16.npz'), _tri('parts_18.npz'), _tri('parts_20.npz'),
        _tri('parts_22.npz'))))
chk('those four partitions really are perfect: 3T = |W|',
    _tri('parts_16.npz') is not None
    and 3 * _tri('parts_16.npz') == _npy_rows(os.path.join(
        KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data', 'lam16_W.npy'))
    and 3 * _tri('parts_18.npz') == _npy_rows(os.path.join(
        KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data', 'lam18_W.npy'))
    and 3 * _tri('parts_20.npz') == _npy_rows(os.path.join(
        KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data', 'lam20_W.npy'))
    and 3 * _tri('parts_22.npz') == _npy_rows(os.path.join(
        KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data', 'lam22_W.npy')),
    'a partition the paper calls perfect does not cover its cap set')
chk('the Kappa sections really are larger than the laminated ones they replaced',
    _capsize(12) > 648 and _capsize(13) > 906,
    'K_12 %s against 648, K_13 %s against 906' % (_capsize(12), _capsize(13)))


# ---- the odd k: the partitions, that they are perfect, and that the TOTAL follows ------
# The paragraph states four old counts, four new ones and a number of spheres.  The last is a
# function of the first eight and of the shipped class sizes, so it is recomputed: a
# sentence's numbers are one calculation, and guarding them separately guards nothing.
def _npy_ints(path):
    """the int64 contents of a .npy, without numpy"""
    import array, ast
    with io.open(path, 'rb') as h:
        raw = h.read()
    i = raw.index(b'{')
    j = raw.index(b'}', i) + 1
    hdr = ast.literal_eval(raw[i:j].decode('latin-1'))
    code = {'<i2': 'h', '<i4': 'i', '<i8': 'q'}.get(hdr['descr'])
    assert code and not hdr['fortran_order'], hdr
    off = raw.index(chr(10).encode(), j) + 1
    a = array.array(code)
    a.frombytes(raw[off:off + a.itemsize * hdr['shape'][0]])
    return list(a)


_CAPD = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
_FAM = []
for _fn in ('class_sizes_32000.npy', 'class_sizes_gpu.npy'):
    _p = os.path.join(_CAPD, _fn)
    if os.path.exists(_p):
        _sz = sorted(_npy_ints(_p), reverse=True)
        _cum = [0]
        for _x in _sz:
            _cum.append(_cum[-1] + _x)
        _FAM.append(_cum)


def _capgain(T, P):
    return max(4 * c[T] + 2 * (c[T + P] - c[T]) for c in _FAM if len(c) > T + P)


# The paragraph runs over several lines, so it is matched against a whitespace-normalised
# copy: a regex with \s+ between every word is a regex nobody can read or repair.
_ODDT = re.sub(r'\s+', ' ', FLAT)
_ODD = re.search(
    r'takes the partitions from \$(' + _num3 + r')\$ triples and \$(' + _num3 +
    r')\$ pairs to \$(' + _num3 + r')\$, from \$(' + _num3 + r')\$ and \$(' + _num3 +
    r')\$ to \$(' + _num3 + r')\$, from \$(' + _num3 + r')\$ and \$(' + _num3 +
    r')\$ to \$(' + _num3 + r')\$, and from \$(' + _num3 + r')\$ and \$(' + _num3 +
    r')\$ to \$(' + _num3 + r')\$, perfect in every case, and moves dimensions '
    r'\$89\$, \$91\$, \$93\$ and \$95\$ by \$(' + _num3 + r')\$ spheres', _ODDT)
def _lamrows(k):
    return _npy_rows(os.path.join(_CAPD, 'lam%d_W.npy' % k))


chk('the paper states the four odd-k partitions, and they are the shipped ones',
    _ODD is not None
    and [num(_ODD.group(i)) for i in (3, 6, 9, 12)]
    == [_tri('parts_%d.npz' % k) for k in (17, 19, 21, 23)],
    'paper %s, shipped %s'
    % ([num(_ODD.group(i)) for i in (3, 6, 9, 12)] if _ODD else None,
       [_tri('parts_%d.npz' % k) for k in (17, 19, 21, 23)]))
chk('all four really are perfect: 3T = |W| and no vector is left in a pair or alone',
    all(3 * _tri('parts_%d.npz' % k) == _lamrows(k) for k in (17, 19, 21, 23)),
    'shortfalls %s'
    % {k: _lamrows(k) // 3 - _tri('parts_%d.npz' % k) for k in (17, 19, 21, 23)})
_want = (sum(_capgain(_tri('parts_%d.npz' % k), 0) for k in (17, 19, 21, 23))
         - sum(_capgain(num(_ODD.group(i)), num(_ODD.group(i + 1)))
               for i in (1, 4, 7, 10))) if (_ODD and _FAM) else None
chk('the spheres the odd-k paragraph claims are what its own old and new counts give',
    _want is not None and num(_ODD.group(13)) == _want,
    'paper %s, recomputed %s' % (num(_ODD.group(13)) if _ODD else None, _want))


# Table 3 needed a set of K(k) poles in each column, and the k = 7 column USED to assume
# one -- the search reached 124 of 126.  The rotated E_7 layer is now all 126, so nothing in
# that table is assumed.  Derive it: every k of 2..7 must have a shipped layer of tau(k), and
# the paper must not still carry the sentence that says otherwise.
def _calib_poles(k):
    """the layer dimension 72+k ships, rotated where there is one and the shell otherwise

    _shipped_poles() insists on a ROTATED layer and exits if there is none; k = 2 and k = 4
    reach tau(k) under the shell rule and have never needed one."""
    import zipfile
    d = os.path.join(KV, 'verifications', 'improved', 'dim73-95-gamma72-caps', 'data')
    for f in ('poles_rot_%d.npz' % k, 'poles_rot_k%d.npz' % k):
        p = os.path.join(d, f)
        if os.path.exists(p):
            with zipfile.ZipFile(p) as z:
                head = z.open('keep.npy').read(128)
            i = head.index(b'(') + 1
            return int(head[i:head.index(b',', i)])
    for f in ('poles_%d.npy' % k, 'poles_k%d.npy' % k):
        p = os.path.join(d, f)
        if os.path.exists(p):
            return 2 * _npy_rows(p)
    return 0


_calib = dict((k, _calib_poles(k)) for k in range(2, 8))
chk('every calibration column exhibits its K(k) poles',
    all(_calib[k] == TAU[k] for k in _calib),
    'short at %s' % [(k, _calib[k], TAU[k]) for k in _calib if _calib[k] != TAU[k]])
# Match the whitespace-collapsed text: a guard that depends on where a line happens to
# wrap fails on a rewrap that changed nothing, and then gets "fixed" by relaxing it.
chk('the paper says so, and no longer says the k=7 column is assumed',
    'is exhibited here in every column' in FLAT
    and 'no column of Table~' in FLAT and 'assumes a pole set' in FLAT
    and 'we have not constructed one' not in FLAT,
    'the calibration sentence does not match the data')

_all = sorted(k for k in _frac if _frac[k][0] == _frac[k][1])
_ceil = sorted(k for k in _frac if _frac[k][1] == TAU[k])
# Read the two k-lists OUT of the paper.  Writing them in here is what made the cap measure
# un-correctable: a literal in a check turns today's value into a requirement.
_said_all = re.findall(r'\d+', (re.search(r'all of \$Z\$ at ([^,.]*(?:, [^,.]*)*?), and reaches',
                                          FLAT) or type('x', (), {'group': lambda s, n: ''})()
                                ).group(1))
_said_ceil = re.findall(r'\d+', (re.search(r'and reaches \$K\(k\)\$ at ([^.]*)\.', FLAT)
                                 or type('x', (), {'group': lambda s, n: ''})()).group(1))
chk('the paper lists the k at which the layer is all of Z, and the list is right',
    [int(x) for x in _said_all] == _all, 'paper %s, data %s' % (_said_all, _all))
chk('the paper lists the k at which the layer reaches K(k), and the list is right',
    [int(x) for x in _said_ceil] == _ceil, 'paper %s, data %s' % (_said_ceil, _ceil))

# ------------------------------------------------------------- prose numbers vs table
print('== Prose numbers against the tables ==')
for d, v in [(25, 197569), (27, 200540), (38, 591612), (39, 756116),
             (70, 1249778250), (71, 2603658750), (96, 12886999232)]:
    chk('prose dim %d = table' % d, tab4[d][1] == v, '%d vs %d' % (tab4[d][1], v))
# Theorem 6.3 and the abstract both state dimension 96.  FLAT collapses whitespace but
# KEEPS the thin spaces, so strip those here.  The theorem writes the bound with no
# space after the relation and the abstract with one, so allow either.
_nothin = TEX.replace(chr(92) + ',', '')
_d96 = re.findall(r'K\(96\)\\ge\s*12886999232', _nothin)
chk('the paper states K(96) >= 12886999232 in words, twice',
    len(_d96) >= 2, 'found %d' % len(_d96))
# A referee meets the ERS abstract before its section 3, and "in dimensions from 32 to 128"
# reads as a claim at every length unless the paper says otherwise.  Require that it does.
# Sun-Wang reconstruct dimensions 38-47 and beat the tabulated bounds in 42-47.  That
# sentence is what protects this paper's only two in-table claims, 38 and 39; without it the
# reader is left to infer it from which dimensions Cohn's table credits to whom.
chk('the paper records the Sun-Wang 42-47 range',
    '$42$--$47$' in TEX and 'only two this paper claims' in FLAT,
    'the reconstruction range is not stated')
chk('the paper confronts the range in the ERS abstract',
    'span of those seven' in FLAT and 'tabulates no value by dimension' in FLAT,
    'the abstract of [ERS] is not addressed')
chk('abstract range 25..96', min(tab4) == 25 and max(tab4) == 96)
chk('abstract count matches Table 4',
    '$%d$ dimensions between $25$ and $96$' % len(tab4) in TEX,
    'the abstract does not say %d' % len(tab4))
chk('title states the range, not a tally',
    'in dimensions 25 through 96' in FLAT and '47 dimensions' not in FLAT)
# 'no published table has an entry' is the claim; d >= 49 is only a proxy for it, and the
# proxy dies quietly the day a published table gains an entry above 49.  Test the claim
# itself against published.py -- which has no imports of its own, so this costs nothing --
# and require the proxy to agree.  A divergence means the table moved, not that the code did.
sys.path.insert(0, os.path.join(KV, 'common'))
try:
    import published as _pub
except Exception as _e:
    chk('the untabulated count can read published.py', False, repr(_e))
    _pub = None
if _pub is not None:
    _untab = sorted(d for d in tab4 if d not in _pub.COHN)
    _proxy = sorted(d for d in tab4 if d >= 49)
    chk('43 of the 47 have no entry in the published table',
        len(_untab) == 43, str(len(_untab)))
    chk('the four that do are 25, 27, 38 and 39',
        sorted(set(tab4) - set(_untab)) == [25, 27, 38, 39],
        str(sorted(set(tab4) - set(_untab))))
    chk('the d >= 49 proxy still agrees with the published table', _untab == _proxy,
        'they differ at %s' % sorted(set(_untab) ^ set(_proxy)))
    # _untab tests membership of COHN, so it cannot see a Nebe-Sloane entry, and the
    # manuscript said twice that no table has an entry above 48.  Dimension 80 is one:
    # Gamma_72 perp E_8.  Both sentences now name it, and both are checked against the
    # dimension list rather than against the wording alone.
    _tab48 = tuple(sorted(d for d in tab4 if d > 48
                          and (d in _pub.NEBE_SLOANE or d in _pub.COHN)))
    chk('exactly one claimed dimension above 48 is in a published table',
        _tab48 == _pub.TABULATED_ABOVE_48 == (80,), str(_tab48))
    chk('and its floor is the Nebe-Sloane entry to the digit',
        tab4[80][0] == _pub.COHN[72] + 240, str(tab4[80][0]))
    _n80 = sum(1 for _s in
               ('ranges in which the published tables have one entry between them, at $n=80$',
                'the tables have a single entry, at $n=80$')
               if _s in FLAT)
    chk('the two no-table sentences both name dimension 80', _n80 == 2,
        '%d of 2 found' % _n80)
    chk('and neither says no table has an entry any more',
        'no published table has an entry' not in FLAT
        and 'no table has an entry' not in FLAT, 'the old wording survives')
    # Section 1.2 has to say where the two summands come from and that the tabular
    # floor is not the only one; a referee meets it long before section 8.2.
    for _s in ("where $K(n-48)$ is Cohn's entry", "$K(n-72)$ is again Cohn's entry",
               'A table entry is not the only thing a bound above dimension $48$',
               'is the tabular floor defined here'):
        chk('section 1.2 states %r' % _s[:42], _s in FLAT, 'missing')
chk('paper says forty, not thirty-eight',
    'thirty-eight' not in TEX.lower(), 'thirty-eight still present')
# Dimension 81 used to fall 101 070 short of dimension 80 and be carried by monotonicity of
# K.  It now has a construction of its own, so the paper states the margin instead, and this
# check reads both numbers out of Table 4 rather than repeating either.
_m81 = re.search(r'exceeds dimension \$80\$ by \$([0-9' + re.escape(BS) + r',]+)\$', TEX)
chk('the paper states a dimension-81 margin', _m81 is not None,
    'the margin sentence has changed shape')
chk('the margin the paper states is the difference Table 4 gives',
    _m81 is not None and num(_m81.group(1)) == tab4[81][1] - tab4[80][1],
    'paper %s, table %d' % (_m81.group(1) if _m81 else None,
                            tab4[81][1] - tab4[80][1]))
chk('Gamma72 point range 1.98e5 .. 6.47e9',
    abs(tab4[25][1] - 1.98e5) < 500 and abs(tab4[95][1] - 6.47e9) < 5e7,
    '%d .. %d' % (tab4[25][1], tab4[95][1]))
chk('dim 47 upper bound 621658419 recorded somewhere', '621\\,658\\,419' in TEX)



# ---------------------------------------------------------- new text introduced by the fix
print('== corrected passages present ==')
chk('calibration no longer pins t=2/3 for every column',
    'so $\\gamma=\\tfrac14$ and $t=\\tfrac23$' not in TEX)
chk('calibration mentions the k=1 level rising to 3/4',
    'raises the level to $t=\\tfrac34$ there' in FLAT)
chk('calibration mentions the 30-degree pole separation',
    'the $30^\\circ$ separation of' in FLAT)
chk('parts row uses a x b', '$12{\\times}3{+}2{\\times}2$' in TEX)
chk('certificate section says maxima round upwards',
    'a certified maximum one above it' in FLAT)
chk('range-stops paragraph names the dimension-96 code',
    '$[96,33,24]$' in TEX and '2^{33}' in TEX)
chk('2^30 percentage right', abs(2 ** 30 / 6218175600 - 0.17) < 0.005,
    '%.4f' % (2 ** 30 / 6218175600))
# the paragraph now ends on the ERS total rather than on a percentage, so the total and
# the factor it quotes are checked against closed/dim96-ers-takeover/ers96.py
chk('ERS(96) total in the paper', '12\\,886\\,999\\,232' in TEX)
chk('ERS(96) factor on Gamma_72', abs(12886999232 / 6218175600 - 2.07) < 0.005,
    '%.4f' % (12886999232 / 6218175600))
# The cap value at k = 24 is computed by cap96.py and copied into table48-96.py, ers96.py
# and run_all.py's label.  Read it from table48-96.py rather than keeping a sixth copy
# here: it went stale once already, when the class family grew and dimension 95 rose
# past it, leaving the paper claiming more at 95 than its own construction gave at 96.
_t4896 = io.open(os.path.join(KV, 'verifications', 'closed', 'dim96-ers-takeover',
                              'table48-96.py'), encoding='utf-8').read()
_cap96 = int(re.search(r'^CAP96 = (\d+)', _t4896, re.M).group(1))
chk('the cap value at k = 24 in the paper matches table48-96.py',
    _grp(_cap96) in TEX, 'repo says %s' % _grp(_cap96))
chk('it is below the ERS total, which is why dimension 96 goes to ERS',
    _cap96 < 12886999232, str(_cap96))
chk('and above the dimension-95 claim, the construction being monotone in n',
    _cap96 > tab4[95][1], '%d vs %d' % (_cap96, tab4[95][1]))
chk('Ozeki width identity 164*94 = 15416', 164 * 94 == 15416)


# ---- why the table stops at 96 ---------------------------------------------------------
# The paragraph after "the layered construction stops winning at 96" states an excess and
# the dimension it occurs in.  Both are downstream of every claim in the paper -- the implied
# floor above 96 is built out of them -- so neither can be typed.  above96.py owns the
# computation and is imported, not re-derived: a second implementation is a second reading,
# and can disagree with the thing it guards.
sys.path.insert(0, os.path.join(KV, 'verifications', 'closed', 'dim96-ers-takeover'))
try:
    from above96 import worst_excess as _a96, TOP as _a96_top
except Exception as _e:                                                # noqa: BLE001
    chk('the above-96 guard can import above96.py', False, repr(_e))
    _a96 = None

if _a96 is not None:
    _exc, _at = _a96()
    _i96 = TEX.find('And the table stops there')
    # strip the LaTeX out of the window so the patterns need no escaping at all
    _w96 = ' '.join(TEX[_i96:_i96 + 1400].split()).replace(chr(92), '').replace('$', '')
    _m96 = re.search('most ([0-9.]+)%, at n=([0-9]+)', _w96)
    chk('the paper says why the table stops at 96', _i96 >= 0 and bool(_m96),
        'the paragraph was not found in the shape this guard reads')
    if _m96:
        chk('the above-96 excess is the one above96.py measures',
            abs(float(_m96.group(1)) - 100 * _exc) < 0.0005,
            'paper %s%%, above96.py %.4f%%' % (_m96.group(1), 100 * _exc))
        chk('the above-96 excess is quoted at the dimension where it occurs',
            int(_m96.group(2)) == _at,
            'paper n = %s, above96.py n = %d' % (_m96.group(2), _at))
    _r96 = re.search('evaluated through n=([0-9]+)', _w96)
    chk('the paper states the range above96.py actually sweeps',
        bool(_r96) and int(_r96.group(1)) == _a96_top,
        'paper %s, above96.py %d' % (_r96.group(1) if _r96 else None, _a96_top))

# ---- section 8.2: the exposure subsection, against ers_exposure.py's own ranking -----------
# Everything section 8.2 says is downstream of a claim, and none of it was checked before:
# 1.47, 8.9 and 4.06 are not integers >= 1000 so unsupported.py never sees them, 45235 passes
# only because a copy sits in common/published.py, and the dimension NAMED is prose.  If a
# claim moves far enough that another dimension becomes the most exposed, the subsection
# argues about the wrong bound with every other check still green.  The ranking is imported,
# never parsed: a parser is a second reading and can disagree with what it guards.
sys.path.insert(0, os.path.join(KV, 'verifications', 'closed', 'dim96-ers-takeover'))
try:
    from ers_exposure import ranked as _ers_ranked
except Exception as _e:
    chk('section 8.2 exposure guard can import the ranking', False, repr(_e))
    _ers_ranked = None

if _ers_ranked is not None:
    _order = _ers_ranked()
    _top, _next = _order[0], _order[1]
    _flat = TEX.replace(chr(92) + ',', '')           # strip LaTeX thin spaces: no backslashes
    _i = _flat.find('most exposed bound in this paper')
    _before, _after = _flat[max(0, _i - 300):_i], _flat[_i:_i + 1200]
    _dims = re.findall(r'[Dd]imension \$(\d+)\$', _before)
    _nxt = re.search(r'next smallest factor is \$([\d.]+)\$, in dimension \$(\d+)\$', _after)
    _mar = re.search(r'margin over the chain is \$([\d.]+)', _after)
    # anchored to A(n, w1, w1), not to being the first '\ge' in the window: a second
    # bound in the same paragraph used to capture this and read 30828 for 45235.
    _REL = r'(?:\\ge|=)'                              # the paper may state either relation
    _thr = re.search(r'A\(%d,%d,%d\)' % (_top[0], _top[4], _top[4]) + _REL + r'([0-9]+)',
                     _after)
    _have = re.search(r'A\(%d,%d,%d\)' % (_top[3][0], _top[4], _top[4]) + _REL + r'([0-9]+)',
                      _after)
    _fac = re.search(r'a factor \$([\d.]+)\$', _after)

    chk('section 8.2 names the most exposed dimension', _i >= 0 and bool(_dims)
        and int(_dims[-1]) == _top[0],
        'paper %s, ers_exposure %s' % (_dims[-1] if _dims else None, _top[0]))
    chk('section 8.2 exposure factor', bool(_fac) and abs(float(_fac.group(1)) - _top[7]) < 0.005,
        'paper %s, ers_exposure %.2f' % (_fac.group(1) if _fac else None, _top[7]))
    chk('section 8.2 level-one threshold', bool(_thr) and int(_thr.group(1)) == _top[6],
        'paper %s, ers_exposure %s' % (_thr.group(1) if _thr else None, _top[6]))
    chk('section 8.2 names the code actually in hand at that weight',
        bool(_have) and int(_have.group(1)) == _top[5],
        'paper %s, ers_exposure %s' % (_have.group(1) if _have else None, _top[5]))
    chk('section 8.2 margin over the chain', bool(_mar)
        and abs(float(_mar.group(1)) / 100 + 1 - _top[1] / float(_top[2])) < 0.0005,
        'paper %s%%, ers_exposure %.4f' % (_mar.group(1) if _mar else None,
                                           _top[1] / float(_top[2])))
    chk('section 8.2 runner-up dimension and factor', bool(_nxt)
        and int(_nxt.group(2)) == _next[0] and abs(float(_nxt.group(1)) - _next[7]) < 0.005,
        'paper %s at %s, ers_exposure %s at %.2f'
        % (_nxt.group(2) if _nxt else None, _nxt.group(1) if _nxt else None,
           _next[0], _next[7]))

    # "every other dimension exceeds 5" -- the third clause of the same sentence, and the
    # one that fails first if a middling claim tightens: it needs the WHOLE tail, not a
    # third data point.
    _rest = re.search(r'every other\s+dimension exceeds \$(\d+)\$', _after)
    _worst = min(_r[7] for _r in _order[2:]) if len(_order) > 2 else None
    chk('section 8.2: every dimension past the runner-up clears the stated floor',
        bool(_rest) and _worst is not None and _worst > float(_rest.group(1)),
        'paper says > %s, the tightest of the remaining %d is %.3f (dimension %s)'
        % (_rest.group(1) if _rest else None, len(_order) - 2,
           _worst if _worst else -1,
           min(_order[2:], key=lambda r: r[7])[0] if len(_order) > 2 else None))


# ---- the dimension-27 direction set: its NAME, derived from the shipped object -------------
# Nothing here guarded a polytope's name.  construction.json always held the A_3 roots and
# every verifier read the file, never the label, so every number stayed right while the prose
# called the set an icosahedron -- past every green check.  A name is a claim about geometry,
# so derive it: read the twelve directions, compute the cosines exactly, and require the paper
# to call them what they are.  The icosahedron is kept in the test as the rival it excludes;
# a check that can only confirm the answer already written is not a check.
_D27 = os.path.join(KV, 'verifications', 'improved', 'dim27-triple-partition',
                    'data', 'construction.json')


def _rmul(a, b):
    """Multiply in Z[sqrt2, sqrt3], coordinates on the basis (1, sqrt2, sqrt3, sqrt6)."""
    return (a[0] * b[0] + 2 * a[1] * b[1] + 3 * a[2] * b[2] + 6 * a[3] * b[3],
            a[0] * b[1] + a[1] * b[0] + 3 * (a[2] * b[3] + a[3] * b[2]),
            a[0] * b[2] + a[2] * b[0] + 2 * (a[1] * b[3] + a[3] * b[1]),
            a[0] * b[3] + a[3] * b[0] + a[1] * b[2] + a[2] * b[1])


def _rip(u, v):
    t = (0, 0, 0, 0)
    for x, y in zip(u, v):
        t = tuple(p + q for p, q in zip(t, _rmul(tuple(x), tuple(y))))
    return t


_dirs = json.load(io.open(_D27, encoding='utf-8'))['cap_directions']
_n2 = _rip(_dirs[0], _dirs[0])
chk('dimension 27: the shipped directions have squared norm 144, rationally',
    _n2 == (144, 0, 0, 0), str(_n2))
_cos, _irr = [], 0
for _i in range(len(_dirs)):
    for _j in range(_i + 1, len(_dirs)):
        _t = _rip(_dirs[_i], _dirs[_j])
        if _t[1] or _t[2] or _t[3]:
            _irr += 1
        _cos.append(F(_t[0], _n2[0]))
chk('dimension 27: every pairwise cosine of the direction set is rational',
    _irr == 0, '%d irrational' % _irr)
_prof = sorted((c, _cos.count(c)) for c in set(_cos))
chk('dimension 27: the twelve directions are a cuboctahedron',
    len(_dirs) == 12 and _prof == [(F(-1), 6), (F(-1, 2), 24), (F(0), 12), (F(1, 2), 24)],
    str(_prof))
# The rival: every non-antipodal icosahedral cosine satisfies 5c^2 = 1, so none of them is
# rational and in particular none is -1/2.  The count below is therefore what excludes the
# icosahedron, and it is counted off the file.  Testing 5c^2 = 1 directly would prove
# nothing: these cosines are Fractions, and no rational satisfies it whatever the data says.
chk('dimension 27: it has the 24 pairs at -1/2 that Corollary 4.3 needs, which an icosahedron cannot',
    _cos.count(F(-1, 2)) == 24, str(_cos.count(F(-1, 2))))

_flat27 = TEX.replace(BS + ',', '')
# Anchor on the theorem's own end marker and walk BACK to the nearest 'Take $Z$ to be':
# searching forward finds the general construction of section 5.3 instead, and the window
# would then span two sections -- wide enough for the right word in the wrong place to
# satisfy it.  (It did, until the perturbation test showed which sentence was being read.)
_ze = _flat27.find('$|W|=K(3)=12$')
_zi = _flat27.rfind('Take $Z$ to be', 0, _ze) if _ze >= 0 else -1
_zsent = _flat27[_zi:_ze] if _zi >= 0 else ''
chk('the paper defines Z where this guard can read it', _zi >= 0 and bool(_zsent))
chk('the paper calls that set a cuboctahedron and the A_3 root system',
    'cuboctahedron' in _zsent and 'A_3' in _zsent,
    ' '.join(_zsent.split())[:70])
chk('the sentence defining Z does not call it an icosahedron', 'cosahed' not in _zsent)
# The word may appear -- it has to, since the paper argues the choice is forced -- but only
# inside the paragraph that rules it out.  Anywhere else it is the old error returning.
_p0 = _flat27.find('The choice of direction set is')
_p1 = _flat27.find(chr(10) + chr(10), _p0 if _p0 >= 0 else 0)
_occ = [_m.start() for _m in re.finditer('cosahed', _flat27)]
chk('every mention of the icosahedron sits inside the paragraph that excludes it',
    _p0 >= 0 and len(_occ) > 0 and all(_p0 < _p < _p1 for _p in _occ),
    '%d mentions, excluding paragraph spans [%d,%d]' % (len(_occ), _p0, _p1))


# ---- every chain tuple the paper names, against the chain finder itself --------------------
# A chain is written (n_0, n_1, ...), and in every dimension here but one n_0 equals n, so
# "(n,15,2)" reads correctly.  Dimension 68 is the exception: its best chain keeps the
# length-64 sign code, so the tuple is (64,16,4,1) and a reader following the convention
# would otherwise go looking for A(68,17).  The tuple is derived here, never parsed from
# a table of expected answers.
if _ers_ranked is not None:
    try:
        from ers_exposure import best_chain as _bc
    except Exception as _e2:
        chk('the chain guard can import best_chain', False, repr(_e2))
        _bc = None
    if _bc is not None:
        for _n, _txt in [(62, '(62,15,2)'), (63, '(63,15,2)'),
                         (68, '(64,16,4,1)'), (96, '(96,24,6,1)')]:
            _val, _ch = _bc(_n)
            _want = '(' + ','.join(str(x) for x in _ch) + ')'
            _inpaper = _want in FLAT.replace(' ', '')
            chk('dimension %d: the paper names the chain the finder returns' % _n,
                _want == _txt and _inpaper,
                'finder %s; in the paper: %s' % (_want, 'yes' if _inpaper else 'NO'))
        # and the two totals the paper prints for those chains
        chk('the chain value quoted for dimension 96 is what the finder computes',
            _bc(96)[0] == 12886999232, str(_bc(96)[0]))
        chk('the chain value quoted for dimension 68 is what the finder computes',
            _bc(68)[0] == 331771800, str(_bc(68)[0]))


# ---- section 1.4's enumeration of the claimed dimensions --------------------------------
# The same shape that had gone stale in CITATION.cff and in the repository README: a prose
# list of the results, which nothing compared with the results.  Read the ranges back.
_i14 = TEX.find('no earlier claim in dimensions')
_sent = TEX[_i14:TEX.find('.', TEX.find('$73$', _i14))] if _i14 >= 0 else ''
_listed = set()
for _a, _b, _c in re.findall(r'\$(\d+)\$\s*-{2,3}\s*\$(\d+)\$|\$(\d+)\$', _sent):
    if _c:
        _listed.add(int(_c))
    else:
        _listed |= set(range(int(_a), int(_b) + 1))
chk('section 1.4 lists exactly the dimensions Table 4 claims',
    _i14 >= 0 and _listed == set(tab4),
    'section 1.4 %s; Table 4 %s (missing %s, extra %s)'
    % (len(_listed), len(tab4), sorted(set(tab4) - _listed) or 'none',
       sorted(_listed - set(tab4)) or 'none'))

# ---- code sizes asserted as equalities --------------------------------------------------
# A(24,6) was an "=" in the dimension-96 proof while the cached Brouwer table brackets it
# 16384-24106, and A(64,16,16) was an "=" twice while the repository has it as ">= 30828".
# Nothing numeric moved -- a larger value only raises an ERS total -- but an equality about
# a quantity nobody knows is still a false statement.  Exact values need a reason; the rest
# must be stated as bounds.
_EXACT = {
    (6, 2): 'A(n,2) = 2^(n-1)',
    (1, 1): 'A(n,1) = 2^n',
    (15, 4): 'single value in the cached Brouwer table, not a bracket',
    (96, 1, 1): 'the weight-one words themselves',
    (68, 4, 4): 'binom(68,3)/4, since 68 = 2 mod 6 admits a Steiner system S(3,4,68)',
}
_eqs = re.findall(r'A\((\d+),(\d+)(?:,(\d+))?\)=', FLAT.replace(BS + ',', ''))
_loose = [tuple(int(x) for x in t if x) for t in _eqs]
_unknown = sorted(set(k for k in _loose if k not in _EXACT))
chk('every code size stated as an equality is one that is known exactly',
    not _unknown,
    'stated with "=" but not exact: %s' % (_unknown or 'none'))
chk('and the exact ones are the ones with a reason on record',
    set(_loose) <= set(_EXACT), str(sorted(set(_loose))))

# ---- claims about Echols's table, against the copy of it in the repository -------------
# "124 new constant-weight codes" and "3323 -> 3324" are facts about an outside paper.  Its
# Table 1 is cached here, so read them off the cache instead of trusting two prose copies.
_ech = os.path.join(KV, 'verifications', 'closed', 'dim32-44-ers-audit', 'data',
                    'echols_table1.txt')
if os.path.exists(_ech):
    _et = io.open(_ech, encoding='utf-8').read()
    _tot = re.search(r'total new bounds:\s*(\d+)', _et)
    chk('the cached Echols table states its own total',
        bool(_tot), 'no total line in echols_table1.txt')
    if _tot:
        chk('the paper reports the same number of new codes',
            ('$%s$ new constant-weight codes' % _tot.group(1)) in FLAT,
            'the cache says %s' % _tot.group(1))
    _398 = re.search(r'A\(39,8,8\):\s*(\d+)\s*->\s*(\d+)', _et)
    chk('the cached table records the A(39,8,8) improvement',
        bool(_398), 'no A(39,8,8) line in echols_table1.txt')
    if _398:
        _lo, _hi = (int(_x) for _x in _398.groups())
        chk('the paper improves A(39,8,8) by the cached amounts',
            ('from $%s$ to $%s$' % (_grp(_lo), _grp(_hi))) in FLAT,
            'the cache says %d -> %d' % (_lo, _hi))
else:
    chk('the Echols table is cached in the repository', False, _ech)

# ---- pointers into other people's documents ---------------------------------------------
# "Equation (3)" was the preprint arXiv:2312.05121's in two files of the repository and the
# journal version's in four; two versions need not number alike, so a pointer has to name the
# one it counts in -- including this comment, which audit check 5f now reads, because the
# paper moved inside the package and its own guards came into scope with it.
# Sun-Wang is a preprint that the repository pins to v3 everywhere else, and the manuscript
# named no version at all.  Neither is a number any recomputation can reach.
for _phrase, _doc in (('equation~(3)', 'arXiv:2312.05121'),
                      ('Table~3 of Sun', 'arXiv:2607.20359v3')):
    _j = TEX.find(_phrase)
    chk('the paper says %r' % _phrase, _j >= 0, 'the pointer is gone')
    chk('and names %s beside it' % _doc,
        _j >= 0 and _doc in TEX[max(0, _j - 300):_j + 300],
        'the pointer does not say which document it numbers')

# The two captions that say which dimensions are absent.  Table 1's used to say 46 and 47
# while Table 4's said 44 to 47, both about the same two tables.  One states it now.
chk('only Table 4 enumerates the absent dimensions',
    'Dimensions $44$ to $47$ are absent:' in FLAT
    and 'Dimensions $44$ to $47$ are absent, for the reasons given with' in FLAT,
    'the captions do not stand in the stated relation')


# ---- the abstract's count of what is NOT lattice-based --------------------------------
# It said "three" before dimension 96 and was corrected by hand; nothing would have caught
# it. tab4's method column is the authority: a bound is one of the lattices' exactly when
# its method names Leech, P_48 or Gamma_72.  The rest -- the constant-weight code at 39 and
# the three ERS chains -- are what "all but" excludes.
_lat = [d for d in tab4 if any(_m in tab4[d][2] for _m in (chr(92) + 'Leech',
                                                           chr(92) + 'Pff',
                                                           chr(92) + 'Gsz'))]
_other = sorted(set(tab4) - set(_lat))
# _NUMWORD, not _WORD: line ~190 already has a capitalised one for the Table 1
# sentence, and a second binding of the same name is a trap even when the
# uses happen not to overlap.
_NUMWORD = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six',
            7: 'seven', 43: 'forty-three'}
_abstract = TEX[TEX.index(chr(92) + 'begin{abstract}'):TEX.index(chr(92) + 'end{abstract}')]
chk('the abstract excludes exactly the bounds no lattice gives',
    ('All but %s' % _NUMWORD.get(len(_other), '?')) in _abstract,
    'derived %d (%s); the abstract does not say "All but %s"'
    % (len(_other), ', '.join(str(d) for d in _other), _NUMWORD.get(len(_other), '?')))
chk('and those are the constant-weight and chain dimensions',
    _other == [39, 62, 63, 96], str(_other))

# "forty-three of the forty-seven" appears twice -- sections 1.2 and 1.5 -- and both
# times counts the dimensions ABOVE 48, which is not the lattice-based set: the two
# partitions of Table 4 happen to have the same size, 43 either way, and reading the
# count off the wrong one passes today and misfires the day a code-theoretic claim is
# added below 49.  The lattice count is the abstract's "All but four", checked above.
_above = sorted(d for d in tab4 if d >= 49)
chk('the manuscript spells the above-48 count the same way',
    '%s of the forty-seven' % _NUMWORD.get(len(_above), '?') in FLAT.lower(),
    'derived %d of %d above dimension 48' % (len(_above), len(tab4)))
chk('and the phrase is used for that count in both places',
    FLAT.lower().count('forty-three of the forty-seven') == 2,
    '%d occurrences' % FLAT.lower().count('forty-three of the forty-seven'))
chk('the two partitions of Table 4 are still distinct sets',
    set(_above) != set(tab4) - set(_other),
    'they coincide, so neither count can be told from the other')

# ERS's seven worked examples: written out in the paper, held as ERS_PUBLISHED's
# keys here.  The paper's whole case for 62, 63 and 96 is that those lengths
# are not in this list, so the list has to be the same in both places.
from published import ERS_PUBLISHED   # noqa: E402
_m = re.search(r'evaluate it only at\s*\$?n\s*=\s*([0-9,\s]+)\$', FLAT)
_ersdims = sorted(int(_x) for _x in re.findall(r'\d+', _m.group(1))) if _m else []
chk('the paper lists exactly the dimensions ERS_PUBLISHED records',
    _ersdims == sorted(ERS_PUBLISHED),
    'paper %s, published.py %s' % (_ersdims, sorted(ERS_PUBLISHED)))
chk('and no chain dimension is among them',
    not ({62, 63, 96} & set(_ersdims)),
    str(sorted({62, 63, 96} & set(_ersdims))))

# Section 6 opens with the same count from the other side -- "Four improvements
# are code-theoretic" -- and said "Three" until now, one section away from
# an abstract that had already been corrected.  Same derived number, so
# check it in the same place.
_sec6 = re.search(r'(\w+) improvements are code-theoretic', TEX)
chk('section 6 counts the code-theoretic improvements the same way',
    bool(_sec6) and _sec6.group(1).lower() == _NUMWORD.get(len(_other), '?'),
    'section 6 says %r; derived %s'
    % (_sec6.group(1) if _sec6 else None, _NUMWORD.get(len(_other), '?')))


# ---- the reference layer -------------------------------------------------------------
# Everything above recomputes a number.  A dangling \\ref carries no number: it renders as
# "??" in the PDF and nothing else notices.  Likewise a .bib entry nobody cites still
# prints a bibliography item, for a work the paper does not use.  Both are read from the
# source, not from .aux or .log, so this holds for the next build and not just the last.
_BIB = io.open(os.path.join(HERE, 'kissing46.bib'), encoding='utf-8').read()
_labels = re.findall(_BS + r'label\{([^}]*)\}', TEX)
_refs = set()
for _m in re.finditer(_BS + r'(?:eq|page|auto|c)?ref\{([^}]*)\}', TEX):
    _refs |= {_k.strip() for _k in _m.group(1).split(',') if _k.strip()}
_dupes = sorted({_x for _x in _labels if _labels.count(_x) > 1})
chk('no label is defined twice', not _dupes, ', '.join(_dupes))
chk('every cross-reference resolves to a label',
    not (_refs - set(_labels)), ', '.join(sorted(_refs - set(_labels))))

_defined = set(re.findall(r'^@[a-zA-Z]+\{([^,]+),', _BIB, re.M))
_cited = set()
for _m in re.finditer(_BS + r'cite[a-zA-Z]*' + r'\s*(?:\[[^]]*\]\s*)*' + r'\{([^}]*)\}', TEX):
    _cited |= {_k.strip() for _k in _m.group(1).split(',') if _k.strip()}
chk('every citation resolves to a bibliography entry',
    not (_cited - _defined), ', '.join(sorted(_cited - _defined)))
chk('every bibliography entry is cited',
    not (_defined - _cited), ', '.join(sorted(_defined - _cited)))
# The .bbl is generated by bibtex and shipped beside the .pdf, because a reviewer who
# builds the paper without running bibtex gets their bibliography from it.  If the .bib
# gained an entry and nobody re-ran bibtex, the two disagree and only this notices.
_bblp = os.path.join(HERE, 'kissing46.bbl')
if os.path.exists(_bblp):
    _bbl = io.open(_bblp, encoding='utf-8').read()
    _inbbl = set(re.findall(_BS + r'bibitem\{([^}]*)\}', _bbl))
    chk('the shipped .bbl lists exactly the entries the .bib defines',
        _inbbl == _defined,
        'only in .bbl: %s; only in .bib: %s'
        % (', '.join(sorted(_inbbl - _defined)) or 'none',
           ', '.join(sorted(_defined - _inbbl)) or 'none'))
else:
    print('   (no kissing46.bbl beside the paper; skipping the stale-bibliography check)')

# The software's title carries a dimension range, and so do the PDF metadata and the .bbl.
# All three are the span of Table 4, which is derived above -- so derive them too.  The .bib
# said "25 through 95" while every other copy said 96, and the stale one is the copy a reader
# of the paper sees.
_span = 'dimensions $%d$ through $%d$' % (min(tab4), max(tab4))
chk('the .bib titles the software by the span of Table 4',
    _span in re.sub(r'\s+', ' ', _BIB),
    'the entry does not say %r' % _span)
# the metadata field itself, not the visible title: the visible one is checked above
# and satisfied this test on its own, so it could not fail for its own reason.
_pdft = re.search(r'pdftitle=\{([^}]*)\}', FLAT)
chk('the PDF metadata title says the same span',
    bool(_pdft) and ('dimensions %d through %d' % (min(tab4), max(tab4)))
    in _pdft.group(1),
    'pdftitle is %r' % (_pdft.group(1) if _pdft else None))
if os.path.exists(_bblp):
    chk('and so does the shipped .bbl',
        ('through %d' % max(tab4)) in re.sub(r'\s+', ' ', _bbl)
        or ('through $%d$' % max(tab4)) in re.sub(r'\s+', ' ', _bbl),
        'the .bbl still carries an older span -- re-run bibtex')

# The repository points at the paper BY SECTION NUMBER, and three files said 7.2 for what
# is 8.2 -- a section was inserted ahead of it and nothing renumbered them.  amsart numbers
# by counting, so count: \section*{...} is unnumbered and this pattern does not match it.
_bd = TEX[TEX.index(BS + 'begin{document}'):]
_ix = _bd.index('label{ss:exposure}')
_sec = len(re.findall(re.escape(BS) + r'section\{', _bd[:_ix]))
_sub = len(re.findall(re.escape(BS) + r'subsection\{',
                      _bd[_bd.rindex(BS + 'section{', 0, _ix):_ix]))
_want = '%d.%d' % (_sec, _sub)
_wrong = []
for _dp, _dn, _fn in os.walk(KV):
    _dn[:] = [_d for _d in _dn if _d not in ('__pycache__', '.git', 'classes', 'state')]
    for _f in _fn:
        if not _f.endswith(('.py', '.md')):
            continue
        _q = os.path.join(_dp, _f)
        try:
            _t = io.open(_q, encoding='utf-8').read()
        except (IOError, UnicodeDecodeError):
            continue
        for _m in re.finditer(r"paper's section (\d+\.\d+)", _t):
            if _m.group(1) != _want:
                _wrong.append((os.path.relpath(_q, KV), _m.group(1)))
chk('the repository points at the right section of the paper',
    not _wrong,
    'section %s is the margin analysis; %s'
    % (_want, '; '.join('%s says %s' % _w for _w in sorted(set(_wrong)))))

# The paper's LENGTH is stated in two READMEs, in two phrasings, and both were stale.
# Only the PDF knows, so read it here and sweep both files from one place.
_pdfp = os.path.join(HERE, 'kissing46.pdf')
if os.path.exists(_pdfp):
    _pages = max(int(_x) for _x in
                 re.findall(b'/Count ([0-9]+)', io.open(_pdfp, 'rb').read()))
    for _where in (os.path.join(HERE, 'README.md'),
                   os.path.join(KV, 'README.md')):
        if not os.path.exists(_where):
            continue
        _t = io.open(_where, encoding='utf-8').read()
        _n = re.findall(r'(\d+) pages,|is a (\d+)-page', _t)
        _n = [int(_a or _b) for _a, _b in _n]
        chk('%s states the paper\'s length'
            % os.path.relpath(_where, os.path.dirname(HERE)).replace(chr(92), '/'),
            bool(_n) and all(_v == _pages for _v in _n),
            'it says %s; the PDF has %d' % (_n or 'nothing', _pages))

# paper/README.md states this script's size in the present tense, and it was 62 checks
# stale.  The number is a property of THIS file and was written in one the file
# never read, so read it.  Counting itself keeps the arithmetic exact.
_rm = os.path.join(HERE, 'README.md')
if os.path.exists(_rm):
    _rmt = io.open(_rm, encoding='utf-8').read()
    _said = re.search(r'`factcheck\.py` \((\d+) checks', _rmt)
    chk('paper/README.md states this script\'s current size',
        bool(_said) and int(_said.group(1)) == OK + 1,
        'README says %s; this run has %d'
        % (_said.group(1) if _said else None, OK + 1))

print()
print('=' * 78)
print('%d checks passed, %d FAILED' % (OK, len(BAD)))
if BAD:
    print()
    for _name, _detail in BAD:
        print('  * %-52s %s' % (_name, _detail))
    print('=' * 78)
    print('*** THE PAPER DISAGREES WITH THE REPOSITORY -- DO NOT SHIP ***')
else:
    print('=' * 78)
    print('ALL CHECKS PASSED')
print('=' * 78)
raise SystemExit(1 if BAD else 0)
