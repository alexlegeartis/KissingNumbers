#!/usr/bin/env python3
"""Rewrite the dimensions 81-95 rows of README.md from the driver's own output.

    python final81-95.py && python scripts/make_readme_table.py

The table had drifted from the driver (it still showed dimension 83 as "carried by dimension
82", and every gain was a few hundred low because the k = 23 axis layer was missing from the
shipped data).  Generating it removes the possibility.
"""
import json, os, re

N = 6218175600

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
rows = json.load(open(os.path.join(PKG, 'rows81-95.json')))
readme = os.path.join(PKG, 'README.md')
text = open(readme, encoding='utf-8', newline='').read()

# Join with the README's own line ending: it is read with newline='' so the surrounding
# lines keep theirs, and a bare newline here would leave a CRLF file mixed.
EOL = '\r\n' if '\r\n' in text else '\n'


def fmt(n):
    return '{:,}'.format(n).replace(',', ' ')


new = []
for r in rows:
    # the source is already the column to the left, so strip it from the scheme -- for
    # EVERY source name, not a list of them: 'kappa' was added on 2026-08-30 and printed
    # twice on the same row until this stopped enumerating.
    scheme = re.sub(r'\((?:lattice|record|record\(Z\[sqrt2\]\)|kappa)\)', '',
                    r['scheme']).replace('  ', ' ').strip()
    hi = '**%s**' % fmt(r['value']) if r['dim'] == 95 else fmt(r['value'])
    w = {'lattice': 'Lambda_k', 'record': 'record', 'record(Z[sqrt2])': 'Z[sqrt2]',
         'kappa': 'K_k'}.get(r['W'], r['W'])
    new.append('| %s | %d | %s | %s | +%s | %s | %s |'
               % ('**95**' if r['dim'] == 95 else r['dim'], r['k'], w, scheme,
                  fmt(r['gain']), hi, fmt(r['previous'])))
new = EOL.join(new)

lines = text.split(EOL)
first = next(i for i, l in enumerate(lines) if l.startswith('| 81 |'))
last = next(i for i, l in enumerate(lines) if l.startswith('| **95** |') or l.startswith('| 95 |'))
out = EOL.join(lines[:first] + new.split(EOL) + lines[last + 1:])

# ---------------------------------------------------------------------------------------
# The FAMILY table and the CROSSOVER table were hand-written and both went stale when the
# greedy family was extended from 31 000 classes to 32 000 on 2026-08-26: the README still
# announced 31 000 classes, sizes 1938-2111 and 63 464 272 lines against the shipped
# 32 000 / 1932-2118 / 65 570 602, and every "fresh" entry of the crossover table was a few
# hundred low.  No claimed VALUE was affected -- final81-95.py reads the .npy files, not the
# README -- but a reader comparing the two would have found them inconsistent.  Both tables
# are generated here from the same files the driver reads.
import numpy as np

DATA = os.path.join(PKG, 'data')
FAM = [('class_sizes_32000.npy', 'automorphism images of one greedy 2106-line class'),
       ('class_sizes_gpu.npy', 'greedy on a 80 000 000-line pool of conflict density '
        '0.004827 (`scripts/gpu/`)')]
sizes, cum, famrows = {}, {}, []
for fn, how in FAM:
    a = np.sort(np.asarray(np.load(os.path.join(DATA, fn))).ravel())[::-1].astype(np.int64)
    sizes[fn], cum[fn] = a, np.concatenate([[0], np.cumsum(a)])
    famrows.append('| `%s` | %s | %d–%d | %s | %s |'
                   % (fn, fmt(len(a)), a.min(), a.max(), fmt(int(a.sum())), how))

lines = out.split(EOL)
h = next(i for i, l in enumerate(lines)
         if l.startswith('| family | classes | sizes | total lines |'))
lines[h + 2:h + 2 + len(famrows)] = famrows

IMG, GPU = FAM[0][0], FAM[1][0]


def img_gain(r):
    """what dimension 72+k reaches on the CPU-regenerable image family alone"""
    c = cum[IMG]
    return (4 * int(c[r['T']]) + 2 * int(c[r['T'] + r['P']] - c[r['T']]) + r['poles'])


cross = []
for r in rows:
    n = r['T'] + r['P']
    if n > min(len(sizes[IMG]), len(sizes[GPU])):
        continue
    a, b = int(cum[IMG][n]), int(cum[GPU][n])
    cross.append('| %d | %s | %s | %s | %+.1f%% |'
                 % (r['dim'], fmt(n), fmt(a), fmt(b), 100.0 * (b - a) / a))
h = next(i for i, l in enumerate(lines)
         if l.startswith('| dim | classes read | images | fresh |'))
e = next(i for i in range(h + 2, len(lines)) if not lines[i].startswith('|'))
lines[h + 2:e] = cross

# What declining the GPU family altogether would cost.  It costs no RECORD: the image
# family is regenerable on a CPU and on its own still clears every published value.
run, margins = 0, []
for r in rows:
    run = max(run, N + img_gain(r))
    margins.append((run - r['previous'], r['dim']))
worst = min(margins)
h = next(i for i, l in enumerate(lines) if l.startswith('`final81-95.py` takes the better'))
import textwrap
_para = ('never mixed, since disjointness ACROSS families is not claimed. Declining the GPU '
         'family entirely costs no RECORD: on the image family alone, which is regenerable '
         'on a CPU, every dimension 81–95 still clears its published value, the '
         'smallest margin being %s at dimension %d.' % (fmt(worst[0]), worst[1]))
_e = next(i for i in range(h + 1, len(lines)) if not lines[i].strip())
lines[h + 1:_e] = textwrap.wrap(_para, 95)   # to the blank line, so re-running is idempotent
# ---------------------------------------------------------------------------------------
# EVERYTHING BELOW DIMENSION 81 was hand-written, and had drifted a pole layer behind the
# driver: dimension 79 read "42 triples + 76 poles / +353 804" against 124 poles and
# +353 852, with the VALUE column beside it already corrected, so the row disagreed with
# itself.  The root-system table was four rows stale and its total, "+100 spheres", was the
# sum of the stale rows.  Generated here from rows73-80.json and the shipped layer files.
import sys

sys.path.insert(0, os.path.join(PKG, '..', '..', '..', 'common'))
from published import TAU                                            # noqa: E402
sys.path.insert(0, HERE)
from axis_rotate import directions, rows as _rows81                  # noqa: E402

SUB = {'1': '\u2081', '2': '\u2082', '3': '\u2083', '4': '\u2084', '5': '\u2085',
       '6': '\u2086', '7': '\u2087', '8': '\u2088'}


def pretty(name):
    """A_3 -> A(subscript 3), so the generated rows read like the hand-written ones did"""
    a, _u, b = name.partition('_')
    return a + ''.join(SUB.get(ch, ch) for ch in b)


def layer_k(k):
    """(shell, rotated or None) for the root system of dimension 72+k"""
    sh = os.path.join(DATA, 'poles_k%d.npy' % k)
    ro = os.path.join(DATA, 'poles_rot_k%d.npz' % k)
    return (2 * len(np.load(sh)) if os.path.exists(sh) else 0,
            len(np.load(ro)['keep']) if os.path.exists(ro) else None)


r7380 = json.load(open(os.path.join(PKG, 'rows73-80.json')))
assert 'table' in r7380, 'run final73-80.py first: rows73-80.json has no table block'
low = []
for r in r7380['table']:
    scheme = r['scheme'].replace('tri +', 'triples +').replace('pair +', 'pairs +')
    scheme = scheme.replace(' 0 pairs +', '').replace('  ', ' ').strip()
    scheme = scheme.replace('1 pairs', '1 pair').replace('1 triples', '1 triple')
    g = '**+%s**' % fmt(r['gain']) if r['dim'] == 80 else '+%s' % fmt(r['gain'])
    low.append('| %d | %d | %s | %s | %s | **%s** | %s |'
               % (r['dim'], r['k'], pretty(r['W']), scheme, g,
                  fmt(r['value']), fmt(r['previous'])))
# `lines` is the live list the family and crossover edits above were made in.  Re-splitting
# `out` here would silently discard them: `out` is only rebuilt at the very end.
f0 = next(i for i, l in enumerate(lines) if l.startswith('| 73 |'))
f1 = next(i for i, l in enumerate(lines) if l.startswith('| 80 |'))
lines[f0:f1 + 1] = low

# ---- the root-system pole table, and the two passages that quote it --------------------
SYSNAME = {3: 'A_3', 4: 'D_4', 5: 'D_5', 6: 'E_6', 7: 'E_7', 8: 'E_8'}
rootrows, gain_shell, at_tau = [], 0, []
for k in range(3, 9):
    sh, ro = layer_k(k)
    use = ro if ro is not None else sh
    gain_shell += use - sh
    if use == TAU[k]:
        at_tau.append(k)
    note = ('already at the ceiling' if ro is None and sh == TAU[k]
            else ('the shell layer is kept' if ro is None else ''))
    if use == TAU[k] and ro is not None:
        note = '**at the ceiling**'
    rootrows.append('| %d | %d | %s | %d | %s | %d | %s |'
                    % (72 + k, k, pretty(SYSNAME[k]), sh,
                       '**%d**' % ro if ro is not None and ro > sh else str(ro or sh),
                       TAU[k], note))
h = next(i for i, l in enumerate(lines) if l.startswith('| dim | k | system | shell |'))
e = next(i for i in range(h + 2, len(lines)) if not lines[i].startswith('|'))
lines[h + 2:e] = rootrows


# ---- the record-vs-lattice cap comparison: the LATTICE column is RECOMPUTED ------------
# It went stale twice without anything noticing: at k = 18 and 20 when scripts/omega_parts.py
# closed those partitions, and at k = 17 and 19 when scripts/oddparts.py closed those.  The
# record column is a measurement from research/dim81-95/recordW2.py and is left as measured;
# the lattice column is just 2T + P of the shipped partition, so it is read from the file.
_h = next(i for i, l in enumerate(lines) if l.startswith('| k | record'))
_e = next(i for i in range(_h + 2, len(lines)) if not lines[i].startswith('|'))
_new = []
for _l in lines[_h + 2:_e]:
    _c = [x.strip() for x in _l.strip().strip('|').split('|')]
    _z = np.load(os.path.join(DATA, 'parts_%d.npz' % int(_c[0])))
    _new.append('| %s | %s | %s | %s | **%d** |'
                % (_c[0], _c[1], _c[2], _c[3], 2 * len(_z['triples']) + len(_z['pairs'])))
lines[_h + 2:_e] = _new


def wrap(para):
    return textwrap.wrap(para, 95)


def andlist(xs):
    xs = [str(x) for x in xs]
    return xs[0] if len(xs) == 1 else ', '.join(xs[:-1]) + ' and ' + xs[-1]


_k3 = layer_k(3)
para = ('\u2014 **+%s spheres** over the shell rule, and the layer reaches \u03c4(k) at '
        'k = %s. At k = 3 the 30\u00b0 caps cover four fifths of S\u00b2, and SAMPLING a '
        'rotation there kept %d against the shell rule\u2019s %d; descent finds all %d, '
        'which is \u03c4(3). A file is written only where the rotation wins, so nothing '
        'regresses. Dimension 74 (A\u2082, six roots, six poles) was already at \u03c4(2).'
        % (fmt(gain_shell), andlist(at_tau), 6, _k3[0], _k3[1]))
h = next(i for i, l in enumerate(lines) if l.startswith('\u2014 **+'))
e = next(i for i in range(h + 1, len(lines)) if not lines[i].strip())
lines[h:e] = wrap(para)

# The second passage quoted the shell layers as though they were still what shipped.
para2 = ('k = 11 has no layer at all: dimension 83 uses the 604-point record realised over '
         '\u2124[\u221a2], whose Gram is irrational, and neither rule is set up for it. '
         '[`scripts/verify_poles.py`](scripts/verify_poles.py) checks every layer exactly, '
         'against the same cap directions each claim uses, and the counts are the `+n '
         'poles` column of the table. Nothing beyond the axis layer is claimed. For k = 1 '
         'no pole direction exists at t = 2/3, so dimension 73 uses the pair scheme at '
         't = 3/4 instead, where \u221a(1\u2212t) = 1/2 makes poles automatic.')
h = next(i for i, l in enumerate(lines) if l.startswith('k = 11 has no layer at all'))
e = next(i for i in range(h + 1, len(lines)) if not lines[i].strip())
lines[h:e] = wrap(para2)

def rot_k(k):
    """the largest non-shell layer at this k, whatever built it

    Four files can hold one: poles_rec_<k>.npz from the record configuration of R^k,
    poles_sqrt2_<k>.npz from the same idea over Z[sqrt2], poles_rot_<k>.npz from a rotated
    copy of the cap directions, and poles_kap_<k>.npz from a rotated copy of the KAPPA cap
    directions at k = 12 and 13.  Only one of the last two exists at any k, because they
    are copies of different cap sets and the dimension uses one cap set.
    """
    got = []
    for fn in ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k, 'poles_kap_%d.npz' % k,
               'poles_rot_%d.npz' % k, 'poles_rot_k%d.npz' % k):
        p = os.path.join(DATA, fn)
        if os.path.exists(p):
            got.append(len(np.load(p)['keep']))
    return max(got) if got else None


def source_kind(k):
    """'record', 'sqrt2' or 'caps': what the layer at this k is a copy OF"""
    best, kind = -1, None
    for fn, nm in (('poles_rec_%d.npz' % k, 'record'), ('poles_sqrt2_%d.npz' % k, 'sqrt2'),
                   ('poles_kap_%d.npz' % k, 'caps'), ('poles_rot_%d.npz' % k, 'caps'),
                   ('poles_rot_k%d.npz' % k, 'caps')):
        p = os.path.join(DATA, fn)
        if os.path.exists(p):
            n = len(np.load(p)['keep'])
            if n > best:
                best, kind = n, nm
    return kind


# ---- and the two k-lists that say WHY a layer falls short ------------------------------
_R81 = _rows81()
short_of_W, loses_to_caps, from_record = [], [], []
for k in [k for k in range(9, 24) if k != 11]:
    W, _c, _m = directions(k, _R81)
    lay = rot_k(k) or 0
    kind = source_kind(k)
    if kind == 'record':
        from_record.append(k)
    elif len(W) < TAU[k]:
        short_of_W.append(k)
    if lay < min(len(W), TAU[k]) or (kind == 'record' and lay < TAU[k]):
        loses_to_caps.append(k)
para3 = ('Where a layer is short of \u03c4(k) it is short for one of two reasons, and both '
         'are measured rather than argued. For k = %s the ceiling \u03c4(k) is larger than '
         'the cap direction set itself, and a rotated copy cannot exceed the set it is a '
         'copy of \u2014 which is why the poles at k = %s are a rotated copy of the RECORD '
         'configuration of R^k instead, and reach within a few dozen of \u03c4(k). For '
         'k = %s the rotation additionally loses a handful of points to the '
         'caps; dimension 38 shows that descent on SO(k) closes exactly this gap, and it is '
         'affordable at the small k and not at k = 23, where one gradient step costs '
         '8.7\u00d710\u2079 inner products.'
         % (andlist(short_of_W), andlist(from_record), andlist(loses_to_caps)))
h = next(i for i, l in enumerate(lines) if l.startswith('Where a layer is short of'))
e = next(i for i in range(h + 1, len(lines)) if not lines[i].strip())
lines[h:e] = wrap(para3)

# ---- the per-k layer table, and the sentence that totals it ---------------------------
# Hand-written, and three states behind: 464 at k = 10 against a shipped 510, 906 at k = 13
# against 904, 2334 at k = 15 against 2340.  Its total was the sum of those.
def shell_k(k):
    for fn in ('poles_%d.npy' % k, 'poles_k%d.npy' % k):
        p = os.path.join(DATA, fn)
        if os.path.exists(p):
            return 2 * len(np.load(p))
    return 0


perk, tot_over_shell, ceil_dims = [], 0, []
for k in range(9, 24):
    sh, ro = shell_k(k), rot_k(k)
    use = ro if ro is not None else sh
    tot_over_shell += use - sh
    if use == TAU[k] and use:
        ceil_dims.append(72 + k)
    frac = ('**at the ceiling**' if use == TAU[k] and use
            else '%.1f%%' % (100.0 * use / TAU[k]))
    perk.append('| %d | %d | %s | %s | %s | %s |'
                % (72 + k, k, fmt(sh), fmt(use), fmt(TAU[k]), frac))
h = next(i for i, l in enumerate(lines)
         if l.startswith('| dim | k | shell layer | rotated layer |'))
e = next(i for i in range(h + 2, len(lines)) if not lines[i].startswith('|'))
lines[h + 2:e] = perk

_ceil = andlist(ceil_dims)
para4 = ('That is **%s more spheres** across these dimensions, and the layer reaches '
         '\u03c4(k) exactly in dimension%s %s. Dimension 86 is the case worth naming: its cap '
         'directions are *the same 1932-point code* dimension 38 uses, vector for vector, so '
         'dimension 38\u2019s rotation \u2014 the one found by descending on SO(14) rather than by'
         % (fmt(tot_over_shell), '' if len(ceil_dims) == 1 else 's', _ceil))
h = next(i for i, l in enumerate(lines) if l.startswith('That is **'))
e = next(i for i in range(h + 1, len(lines)) if not lines[i].strip())
lines[h:e] = wrap(para4)

out = EOL.join(lines)
open(readme, 'w', encoding='utf-8', newline='').write(out)
print("README.md rows 81-95 regenerated from rows81-95.json (%d rows)" % len(rows))
print("README.md rows 73-80 regenerated from rows73-80.json (%d rows)" % len(low))
print("README.md root-system table regenerated: +%d over the shell rule, tau(k) at k = %s"
      % (gain_shell, andlist(at_tau)))
print("README.md per-k layer table regenerated: +%d over the shell rule in 81-95, at the "
      "ceiling in %s" % (tot_over_shell, andlist(ceil_dims)))
print("README.md family and crossover tables regenerated from the .npy files")
