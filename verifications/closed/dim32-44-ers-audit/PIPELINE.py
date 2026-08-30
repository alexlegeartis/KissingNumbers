#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
PIPELINE.py -- dimensions 32..44 via the Edel-Rains-Sloane construction,
               fed with the best CURRENTLY PUBLISHED code tables.

Run:  python PIPELINE.py            (pure integer arithmetic, no third-party packages)

=====================================================================================
WHAT THIS DOES
=====================================================================================

The Edel-Rains-Sloane (ERS) construction, "On kissing numbers in dimensions 32 to
128", Electron. J. Combin. 5 (1998) #R22:

    tau_n  >=  N(n)  =  sum_{nu=0..mu}  A(n, n_nu, n_nu) * A(n_nu, ceil(n_nu/4))

taken over a chain of support sizes

    n >= n_0 >= 4*n_1 >= 16*n_2 >= 64*n_3 >= ...          (i.e. n_nu >= 4*n_{nu+1})

Level nu holds the vectors of shape (+-a_nu)^{n_nu} 0^{n-n_nu} with a_nu =
sqrt(n_0/n_nu): every centre has squared norm n_0, the supports come from a binary
constant-weight code C(n, n_nu, n_nu) and the signs from a binary code
C(n_nu, ceil(n_nu/4)).  (Constant-weight distances are even, so "A(n,w,w)" here
always means A(n, 2*ceil(w/2), w).)

Geometry, re-derived (all centres have squared norm n_0, so the 60-degree condition
is <x,y> <= n_0/2):
  * same level, same support : <x,y> = n_nu - 2*d_H(signs) <= n_nu/2
                               <=> d_H >= n_nu/4  <=> d >= ceil(n_nu/4).
  * same level, different supports of overlap m : <x,y> <= m, need m <= n_nu/2,
                               i.e. constant-weight distance 2*ceil(n_nu/2).
  * levels u < v (n_u > n_v) : <x,y> <= a_u a_v n_v = n_0 sqrt(n_v/n_u) <= n_0/2
                               <=> n_u >= 4 n_v.        <-- the chain condition

The formula is confirmed three ways:
  (i)  it reproduces the two large published ERS totals to the digit
       (n=64 and n=80, checked below);
  (ii) it reproduces every entry of H. Cohn's kissing table in dimensions 32..44
       except n=43 (which is held by a different construction) -- checked below;
  (iii) W. Echols, arXiv:2608.13906 (Aug 2026), states the same formula verbatim
       ("tau_n >= sum_{nu=0}^{mu} A(n,n_nu,n_nu) A(n_nu, ceil(n_nu/4))") and uses
       the fixed chain (32,8,2).

=====================================================================================
DATA FETCHED (all files live in ./data/, fetched 2026-08-21)
=====================================================================================

  Andw.html            https://aeb.win.tue.nl/codes/Andw.html
                       A. E. Brouwer, "Bounds for binary constant weight codes".
                       Byte-identical to the copy this project fetched 2026-08-19.
                       Carries W. Echols' new codes (emails 2026-07-20/-25/-29/-08-08,
                       arXiv:2608.13906) -- verified against the paper's Table 1.
                       Main tables cover n <= 32; separate sub-tables cover
                       A(n,4,3), A(n,4,4), A(n,6,5), A(n,6,6), A(n,8,5..8) to n = 64.

  LRS_tableand.html    S. Litsyn, E. M. Rains, N. J. A. Sloane, "Table of Nonlinear
                       Binary Codes" (lower bounds on A(n,d)), snapshot "last updated
                       November 24, 1999".  The live copy at
                       http://www.eng.tau.ac.il/~litsyn/tableand/ now returns HTTP 403.
                       This is still the reference table for unrestricted binary codes
                       of length > 28: Brouwer's own A(n,d) table
                       (https://aeb.win.tue.nl/codes/binary-1.html) stops at n = 28,
                       and Erik Agrell's table (https://codes.se/bounds/unr.html) stops
                       at n = 28 and states "the lower bounds ... were copied from
                       [LRS] in January 2001 and are not regularly updated".

  bklc_*.html          M. Grassl, https://codetables.de/  "Bounds on the minimum
                       distance of linear codes over GF(2)" -- best known LINEAR codes,
                       fetched live 2026-08-21.  Gives A(n,d) >= 2^k lower bounds.

  binary-1.html        Brouwer's A(n,d) table for n <= 28.
  echols.html          arXiv:2608.13906v1 (HTML), used to audit Andw.html.
  cohn_kissing.html    https://cohn.mit.edu/kissing-numbers/ , fetched live 2026-08-21
                       -- the record column this run is measured against.

Only LOWER bounds are ever used as code sizes.  Ranges "a-b" in the tables are read
as the lower end a.  Johnson upper bounds are computed separately and used only to
rule chains OUT.

=====================================================================================
RESULT OF THIS RUN (2026-08-21)
=====================================================================================

NO NEW RECORD.  Maximised over ALL admissible chains with every currently published
lower bound, the ERS construction reproduces the best known kissing number in every
dimension 32..44 and improves none of them:

   n    best chain     N(n)        Cohn table    vs Cohn   best known    verdict
  32    (32,8,2)         346432        345408      +1024       346432    no gain
  33    (32,8,2)         362048        360640      +1408       362048    no gain
  34    (32,8,2)         381124        380868       +256       381124    no gain
  35    (32,8,2)         409548        409548         +0       409548    no gain
  36    (32,8,2)         484568        484568         +0       484568    no gain
  37    (32,8,2)         496232        494312      +1920       496232    no gain
  38    (38,8,2)         570236        566652      +3584       591612    no gain
  39    (39,8,2)         756116        755988       +128       756116    no gain
  40    (40,8,2)        1064368       1064368         +0      1064368    no gain
  41    (40,8,2)        1170384       1170384         +0      1170384    no gain
  42    (40,8,2)        1250676       1250676         +0      1250676    no gain
  43    (43,8,2)        1745692       2060399    -314707      2060399    no gain
  44    (44,8,2)        2948552       2948552         +0      2948552    no gain

The "+" column against Cohn's table is not new: dims 32,33,34,37 are W. Echols'
already-published gain (arXiv:2608.13906), and dims 38,39 are this project's own
earlier gain (../../improved/dim39-ers-constant-weight, and for dim 38 the much
larger ../../improved/dim38-leech-large-codimension, which supersedes it).  Cohn's
table is simply stale in those six dimensions.

WHY THERE IS NOTHING LEFT TO FEED IN

  * A(n,8,8), n = 32..44 -- Brouwer's live page is byte-identical to the copy taken
    on 2026-08-19 and already carries every one of Echols' 124 new codes (his emails
    of 2026-07-20/-25/-29 and 2026-08-08 are cited in the page source, and the page's
    values were audited here against Table 1 of arXiv:2608.13906v1).  Echols improved
    A(n,8,8) only for n in {27,30,31,32,33,34,37,38,39} -- nothing for n = 35,36 or
    40..44, whose values (2157, 2742, 3683, 4510, 5136, 5418, 6622) are exactly the
    ones already used by Cohn's 35,36,40,41,42,44 entries.

  * A(n_0, ceil(n_0/4)) -- the top-level sign codes.  Litsyn-Rains-Sloane (1999) and
    Grassl's live codetables.de agree on every value that matters:
        A(32, 8) = 131072 = 2^17   LRS group code [chen89] / linear [32,17,8]
        A(36, 9) =  98304 = 3*2^15 LRS shortened QR [sloa72a]        (< 131072)
        A(38,10) = 180224 = 11*2^14 Zinov'ev-Litsyn Y1 [zino84a]
        A(39,10) = 327680 =  5*2^16 Zinov'ev-Litsyn Y1 [zino84a]
        A(40,10) = 589824 =  9*2^16 Zinov'ev-Litsyn Y1 [zino84a]
        A(41,11) = 262144 = 2^18   shortened QR[47,24,11] / linear [41,18,11]
        A(42,11) = 524288 = 2^19   shortened QR[47,24,11] / linear [42,19,11]
        A(43,11) =1048576 = 2^20   shortened QR[47,24,11] / linear [43,20,11]
        A(44,11) =2097152 = 2^21   shortened QR[47,24,11] / linear [44,21,11]
    Both A(41,11) and A(42,11) sit BELOW A(40,10) = 589824, which is why the optimal
    chain for n = 41 and n = 42 has n_0 = 40, not n_0 = n.  Neither Brouwer's own
    A(n,d) table nor Erik Agrell's goes past n = 28, and Agrell's states that its
    lower bounds were copied from LRS in January 2001 and are not maintained -- so
    LRS + codetables.de really is the current state of the art here.
    codetables.de moreover gives d(42,20) = d(43,21) = d(44,22) = 10 EXACTLY (not a
    range), so no linear code can push A(42,11), A(43,11) or A(44,11) any higher;
    and the open cells [36,17] (8-9) and [37,17] (9-10) would at best TIE 131072.

  * middle levels of weight 9, 10, 11 (admissible once n_0 >= 36, 40, 44).  For w = 9
    the Johnson bound rules them out outright in every dimension.  For w = 10 and 11
    no lower bound on A(n,10,10) or A(n,12,11) is published anywhere for n > 32:
    Brouwer's largest tabulated values are A(32,10,10) >= 500 and A(32,12,11) >= 186,
    while beating the w = 8 level at n = 44 would need A(44,10,10) >= 11773 or
    A(44,12,11) >= 5887 -- 24x and 32x the largest published value, against Johnson
    ceilings of 29682 and 13964.

  * chains with n_0 < 32 lose the weight-8 level entirely and cap the top term at
    A(31,8) = 65536; the enumeration covers them and none comes close.

WHAT WOULD MAKE A RECORD (see section [7] of the output)

  Every dimension except 38 and 43 is exactly at the ERS value, so ANY improvement,
  even by one codeword, is a new record:
     +1 in A(n,8,8) gives +128 in tau_n for n = 32..37, 39..42, 44
     +1 in A(32,8) gives +1 in tau_n for n = 32..37 simultaneously
     +1 in A(40,10) gives +1 in tau_n for n = 40, 41 and 42 simultaneously
     +1 in A(44,11) gives +1 in tau_44
  Dimension 38 needs +21377 (A(38,8,8) >= 3193 or A(38,10) >= 201601) to beat this
  project's Leech-route 591612; dimension 43 needs +314708 (A(43,8,8) >= 7877 or
  A(43,11) >= 1363284) to beat the 2060399 held by Sun-Wang / the P_48 cross-section.
"""

import os, re, sys, html, json
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

def load(fn, enc="utf-8"):
    return open(os.path.join(DATA, fn), encoding=enc, errors="replace").read()

FAIL = []
def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        FAIL.append(msg)

# =====================================================================================
# 1.  A(n,d) -- unrestricted binary codes.  LOWER BOUNDS ONLY.
# =====================================================================================

BIN = {}          # (n,d) -> [value, provenance]
def putbin(n, d, v, prov):
    if v is None or v < 1:
        return
    if (n, d) not in BIN or v > BIN[(n, d)][0]:
        BIN[(n, d)] = [v, prov]

NMAX = 60

# ---- trivial / exact facts -----------------------------------------------------
for n in range(1, NMAX + 1):
    putbin(n, 1, 1 << n, "trivial A(n,1)=2^n")
    putbin(n, 2, 1 << (n - 1), "trivial A(n,2)=2^(n-1)")
    for d in range(1, NMAX + 1):
        if d > n:
            putbin(n, d, 1, "trivial d>n")
        elif 3 * d > 2 * n:
            putbin(n, d, 2, "trivial 3d>2n")

# ---- Brouwer's table for n <= 28 (even d) --------------------------------------
_h = load("binary-1.html")
_m = re.search(r"<TABLE BORDER=\"1\">(.*?)</TABLE>", _h, re.S | re.I)
_rows = re.findall(r"<TR>(.*?)(?=<TR>|</TABLE>)", _m.group(1), re.S | re.I)
_ds = None
for _r in _rows:
    _c = [html.unescape(re.sub(r"<[^>]+>", "", x)).strip() for x in
          re.findall(r"<TD>(.*?)</TD>", _r, re.S | re.I)]
    _c = [x for x in _c if x != ""]
    if not _c:
        continue
    if _c[0].startswith("d="):
        _ds = [int(x.split("=")[1]) for x in _c]
        continue
    if _ds and _c[0].isdigit():
        _n = int(_c[0])
        for _d, _v in zip(_ds, _c[1:]):
            _mm = re.match(r"^(\d+)", _v.replace(" ", ""))
            if _mm:
                putbin(_n, _d, int(_mm.group(1)), "Brouwer binary-1.html")

# ---- Litsyn-Rains-Sloane (odd d, n up to 512) ----------------------------------
_s = load("LRS_tableand.html", "latin-1")
def _lrs_table(d):
    i = _s.find('NAME="dist%d"' % d)
    j = _s.find('NAME="dist%d"' % (d + 2))
    seg = _s[i:j if j > 0 else len(_s)]
    rows = []
    for m in re.finditer(r"<tr>(.*?)(?=<tr>|</table>|$)", seg, re.S | re.I):
        c = re.findall(r"<t[dh]>(.*?)</t[dh]>", m.group(1), re.S | re.I)
        if len(c) < 6:
            continue
        c = [html.unescape(re.sub(r"<[^>]+>", "", x)).strip() for x in c]
        if not c[0].isdigit():
            continue
        rows.append((int(c[0]), int(c[2]), int(c[3]), c[4], c[5]))
    return rows
LRS_RAW = {}
for _d in (3, 5, 7, 9, 11, 13, 15, 17):
    for (_n, _k, _N, _ty, _ref) in _lrs_table(_d):
        LRS_RAW[(_n, _d)] = (_N << _k, _ty, _ref)
        putbin(_n, _d, _N << _k, "LRS 1999 listed %s %s" % (_ty, _ref))
# shortening: a listed (n0,d) entry of size N*2^k yields A(n0-s,d) >= N*2^(k-s)
for (_n0, _d), (_v, _ty, _ref) in list(LRS_RAW.items()):
    _x = _v
    for _s2 in range(1, 16):
        _x //= 2
        if _x < 1:
            break
        putbin(_n0 - _s2, _d, _x, "LRS 1999 shortened from n=%d (%s %s)" % (_n0, _ty, _ref))

# ---- codetables.de best known LINEAR codes -------------------------------------
def _bklc(fn):
    h = load(fn, "latin-1")
    m = re.search(r"<TABLE BORDER=1 CELLPADDING=2.*?</TABLE>", h, re.S | re.I)
    rows = re.findall(r"<TR>(.*?)</TR>", m.group(0), re.S | re.I)
    hdr = [re.sub(r"<[^>]+>", "", c).strip()
           for c in re.findall(r"<TD[^>]*>(.*?)</TD>", rows[0], re.S | re.I)]
    ks = [int(x) for x in hdr[1:]]
    out = {}
    for r in rows[1:]:
        cells = [re.sub(r"<[^>]+>", "", c).strip().replace("&nbsp;", "")
                 for c in re.findall(r"<TD[^>]*>(.*?)</TD>", r, re.S | re.I)]
        if not cells or not cells[0].isdigit():
            continue
        n = int(cells[0])
        for k, c in zip(ks, cells[1:]):
            mm = re.match(r"^(\d+)(?:\s*[-–]\s*(\d+))?$", c.strip())
            if mm:
                out[(n, k)] = int(mm.group(1))     # LOWER end: an actual code exists
    return out
BKLC = {}
for _fn in ("bklc_1_50_lowk.html", "bklc_1_50_hik.html", "bklc_30_50.html"):
    if os.path.exists(os.path.join(DATA, _fn)):
        for _key, _v in _bklc(_fn).items():
            BKLC[_key] = max(BKLC.get(_key, 0), _v)
for (_n, _k), _d in BKLC.items():
    for _dd in range(1, _d + 1):
        putbin(_n, _dd, 1 << _k, "codetables.de linear [%d,%d,%d]" % (_n, _k, _d))

# ---- closure: monotone in n and d, and A(n,2e-1) = A(n+1,2e) --------------------
for _it in range(6):
    for _n in range(1, NMAX + 1):
        for _d in range(1, NMAX + 1):
            if (_n, _d) in BIN:
                v, p = BIN[(_n, _d)]
                putbin(_n + 1, _d, v, "monotone in n from " + p)         # A(n+1,d)>=A(n,d)
                putbin(_n, _d - 1, v, "monotone in d from " + p)         # A(n,d-1)>=A(n,d)
                if _d % 2 == 1:
                    putbin(_n + 1, _d + 1, v, "A(n,2e-1)=A(n+1,2e) from " + p)
                else:
                    putbin(_n - 1, _d - 1, v, "A(n,2e)=A(n-1,2e-1) from " + p)

def Abin(n, d):
    return BIN[(n, d)][0]
def Abin_src(n, d):
    return BIN[(n, d)][1]

# =====================================================================================
# 2.  A(n,d,w) -- binary constant weight codes.  LOWER BOUNDS ONLY.
# =====================================================================================

CW = {}           # (n,d,w) -> [value, provenance]
def putcw(n, d, w, v, prov):
    if v is None or v < 1:
        return
    if (n, d, w) not in CW or v > CW[(n, d, w)][0]:
        CW[(n, d, w)] = [v, prov]

_A = load("Andw.html")

# ---- main tables (n <= 32), one per distance -----------------------------------
for _chunk in re.split(r'<h1><a name="d(\d+)">', _A)[1:]:
    pass
_parts = re.split(r'<h1><a name="d(\d+)">', _A)
for _i in range(1, len(_parts), 2):
    _d = int(_parts[_i]); _body = _parts[_i + 1]
    _m = re.search(r"<table[^>]*>(.*?)</table>", _body, re.S)
    if not _m:
        continue
    _hdr = None
    for _r in re.findall(r"<tr>(.*?)</tr>", _m.group(1), re.S):
        _cells = re.findall(r"<t([hd])[^>]*>(.*?)</t[hd]>", _r, re.S)
        if not _cells:
            continue
        if _cells[0][0] == "h" and _hdr is None and "n\\w" in _cells[0][1]:
            _hdr = [int(re.sub(r"<[^>]+>", "", c[1]).strip()) for c in _cells[1:]]
            continue
        if _cells[0][0] != "h" or _hdr is None:
            continue
        try:
            _n = int(re.sub(r"<[^>]+>", "", _cells[0][1]).strip())
        except ValueError:
            continue
        for _j, (_t, _c) in enumerate(_cells[1:]):
            if _j >= len(_hdr):
                continue
            _w = _hdr[_j]
            _txt = re.sub(r"<sup>.*?</sup>", "", _c, flags=re.S)
            _txt = html.unescape(re.sub(r"<[^>]+>", "", _txt)).strip().rstrip(".")
            _mm = re.match(r"^\s*(\d+)\s*(?:[-–].*)?$", _txt)
            if _mm:
                putcw(_n, _d, _w, int(_mm.group(1)), "Brouwer Andw d=%d main table" % _d)

# ---- the large-n sub-tables ("Values of A(n,d,w)" / "Lower bounds for ...") -----
for _sec in re.split(r"<h3>", _A):
    _title = re.sub(r"<[^>]+>", "", _sec.split("</h3>")[0]).strip()
    _tm = re.search(r"A\(n,(\d+),(\d+)\)", _title)
    if not _tm:
        continue
    _d, _w = int(_tm.group(1)), int(_tm.group(2))
    _also = re.search(r"A\(n,(\d+),(\d+)\)\s+and\s+A\(n,(\d+),(\d+)\)", _title)
    _tab = re.search(r"<table.*?</table>", _sec, re.S | re.I)
    if not _tab:
        continue
    _ns = []
    _wcol = None
    for _r in re.findall(r"<tr>(.*?)</tr>", _tab.group(0), re.S | re.I):
        _cells = [re.sub(r"<sup>.*?</sup>", "", x, flags=re.S)
                  for x in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", _r, re.S | re.I)]
        _cells = [html.unescape(re.sub(r"<[^>]+>", "", x)).strip() for x in _cells]
        _cells = [x for x in _cells if x != "" and x != "n"]
        if not _cells:
            continue
        _isn = (all(re.fullmatch(r"\d+", x) for x in _cells) and len(_cells) >= 4
                and int(_cells[-1]) - int(_cells[0]) == len(_cells) - 1)
        if _isn:
            _ns = [int(x) for x in _cells]
            _wcol = None
            continue
        if not _ns:
            continue
        if _cells[0].startswith("w="):          # the A(n,4,3)/A(n,4,4) combined table
            _wcol = int(_cells[0].split("=")[1]); _vals = _cells[1:]
        else:
            _vals = _cells
        _ww = _wcol if _wcol is not None else _w
        for _n, _v in zip(_ns, _vals):
            _v = _v.strip().rstrip(".")
            _mm = re.match(r"^(\d+)(?:\s*[-–].*)?$", _v)
            if _mm:
                putcw(_n, _d, _ww, int(_mm.group(1)),
                      "Brouwer Andw sub-table '%s'" % _title)
        if _wcol is None:
            _ns = []

# ---- trivial constant-weight facts ---------------------------------------------
for _n in range(1, NMAX + 1):
    putcw(_n, 2, 1, _n, "trivial A(n,2,1)=n")
    putcw(_n, 2, 2, comb(_n, 2), "trivial A(n,2,2)=C(n,2)")
    for _w in range(1, 20):
        putcw(_n, 2 * _w, _w, _n // _w, "trivial: disjoint supports")

# ---- closure: A(n,d,w) is non-decreasing in n ----------------------------------
for _it in range(NMAX):
    for (_n, _d, _w) in list(CW.keys()):
        v, p = CW[(_n, _d, _w)]
        putcw(_n + 1, _d, _w, v, "monotone in n from " + p)

def Acw(n, w):
    """Lower bound on the level-nu support code: A(n, 2*ceil(w/2), w)."""
    d = 2 * ((w + 1) // 2)
    return CW.get((n, d, w), [1, "no published code; A>=1 trivially"])[0]
def Acw_src(n, w):
    d = 2 * ((w + 1) // 2)
    return CW.get((n, d, w), [1, "no published code; A>=1 trivially"])[1]

def johnson_cw(n, d, w):
    """Johnson's nested-floor UPPER bound on A(n,d,w), d even, delta = d/2."""
    if w > n:
        return 0
    delta = d // 2
    if w < delta:
        return 1
    b = (n - w + delta) // delta
    for i in range(1, w - delta + 1):
        b = ((n - w + delta + i) * b) // (delta + i)
    return b

# =====================================================================================
# 3.  THE ERS COUNT, AND THE EXACT MAXIMISATION OVER ALL ADMISSIBLE CHAINS
# =====================================================================================

def sign_dist(w):
    """Minimum distance the level-nu sign code must have."""
    return -((-w) // 4)                       # ceil(w/4)

def term(n, w, cw=None):
    """A(n,w,w) * A(w, ceil(w/4))  -- the contribution of a level of weight w."""
    return (Acw(n, w) if cw is None else cw) * Abin(w, sign_dist(w))

def chains(n, wmax=None, top=True):
    """All admissible chains n >= n_0 >= 4 n_1 >= 16 n_2 >= ... as tuples."""
    hi = n if wmax is None else wmax
    out = []
    for w in range(1, hi + 1):
        out.append((w,))
        for tail in chains(n, w // 4, False):
            out.append((w,) + tail)
    return out

def N_of_chain(n, ch):
    return sum(term(n, w) for w in ch)

def best_chain(n):
    b, bc = -1, None
    for ch in chains(n):
        v = N_of_chain(n, ch)
        if v > b:
            b, bc = v, ch
    return b, bc

# =====================================================================================
# 4.  RUN
# =====================================================================================

def main():
    vol7 = sum(comb(36, i) for i in range(8))
    print("=" * 88)
    print("DIMENSIONS 32-44 VIA EDEL-RAINS-SLOANE, WITH THE CURRENT CODE TABLES")
    print("=" * 88)

    # ---------------------------------------------------------------- sanity checks
    print("\n[1] The formula reproduces the published ERS totals (their own inputs):")
    n64 = (1 << 28) + 30828 * Abin(16, 4) + 10416 * Abin(4, 1) + 64 * Abin(1, 1)
    n80 = (1 << 30) + 143780 * Abin(16, 4) + 20540 * Abin(4, 1) + 80 * Abin(1, 1)
    check(Abin(16, 4) == 2048, "A(16,4) = 2048")
    check(n64 == 331737984, "n=64, chain (64,16,4,1): %d == 331737984" % n64)
    check(n80 == 1368532064, "n=80, chain (80,16,4,1): %d == 1368532064" % n80)

    print("\n[2] Key code sizes now in hand (LOWER bounds, with source):")
    for n0 in range(32, 45):
        d = sign_dist(n0)
        print("      A(%2d,%2d) >= %9d    %s" % (n0, d, Abin(n0, d), Abin_src(n0, d)))
    print()
    for n in range(32, 45):
        print("      A(%2d,8,8) >= %6d    %s" % (n, Acw(n, 8), Acw_src(n, 8)))

    print("\n[3] The three level terms, per dimension (exact integers):")
    print("      %3s %12s %8s %12s %10s" % ("n", "max top", "n_0*", "128*A(n,8,8)", "4*C(n,2)"))
    for n in range(32, 45):
        tops = [(Abin(n0, sign_dist(n0)), n0) for n0 in range(32, n + 1)]
        t, n0 = max(tops)
        print("      %3d %12d %8d %12d %10d" % (n, t, n0, 128 * Acw(n, 8), 4 * comb(n, 2)))

    # ---------------------------------------------------------------- chain search
    print("\n[4] Exhaustive maximisation over ALL admissible chains, n = 32..44")
    print("    (%d chains are enumerated at n = 44)" % len(chains(44)))
    RECORD = {32: 345408, 33: 360640, 34: 380868, 35: 409548, 36: 484568,
              37: 494312, 38: 566652, 39: 755988, 40: 1064368, 41: 1170384,
              42: 1250676, 43: 2060399, 44: 2948552}          # cohn.mit.edu, 2026-08-21
    PROJECT = dict(RECORD)                                    # this project's own floors
    PROJECT.update({32: 346432, 33: 362048, 34: 381124, 37: 496232,
                    38: 591612, 39: 756116})
    # Two of those are this repository's own claims, and a hand-copied claim goes stale the
    # first time the claim moves -- dimension 38 did, on 2026-08-27.  Read RESULTS.md and
    # insist they still agree.
    _res = os.path.join(HERE, "..", "..", "..", "RESULTS.md")
    if os.path.exists(_res):
        for _line in open(_res, encoding="utf-8"):
            _m = re.match(r"\| (\d+) \| [^|]+ \| \*\*([\d\s  ]+)\*\* \|", _line)
            if _m and int(_m.group(1)) in PROJECT:
                _d, _v = int(_m.group(1)), int(re.sub(r"\D", "", _m.group(2)))
                if PROJECT[_d] != _v:
                    raise SystemExit("PIPELINE.py has %d for dimension %d; RESULTS.md says "
                                     "%d -- update the PROJECT dict" % (PROJECT[_d], _d, _v))
    print()
    print("      %3s  %-12s %12s %12s %9s %12s %-12s" %
          ("n", "best chain", "N(n)", "Cohn table", "vs Cohn", "best known", "verdict"))
    print("      " + "-" * 88)
    results = {}
    for n in range(32, 45):
        v, ch = best_chain(n)
        results[n] = (v, ch)
        verdict = ("NEW RECORD +%d" % (v - PROJECT[n])) if v > PROJECT[n] else "no gain"
        print("      %3d  %-12s %12d %12d %+9d %12d %-12s" %
              (n, "(" + ",".join(map(str, ch)) + ")", v, RECORD[n],
               v - RECORD[n], PROJECT[n], verdict))
    print("      " + "-" * 88)
    print("      'vs Cohn' > 0 in dims 32,33,34,37 is Echols' already-published gain")
    print("      (arXiv:2608.13906), and in dims 38,39 this project's own earlier gain;")
    print("      none of it is new here.  'best known' is the max of Cohn's column and")
    print("      the values this project already holds.")

    print("\n[5] Re-derivation of each optimum, term by term (exact):")
    for n in range(32, 45):
        v, ch = results[n]
        parts = []
        tot = 0
        for w in ch:
            t = term(n, w)
            tot += t
            parts.append("%d*%d" % (Acw(n, w), Abin(w, sign_dist(w))))
        print("      n=%2d  %s = %s = %d" % (n, " + ".join(parts),
              " + ".join(str(term(n, w)) for w in ch), tot))
        check(tot == v, "n=%d total re-derives" % n)

    # ------------------------------------------------- why w = 9,10,11 cannot help
    print("\n[6] Could a middle level of weight 9, 10 or 11 beat weight 8?")
    print("    (weight w needs n_0 >= 4w, so w=9 needs n>=36, w=10 needs n>=40,")
    print("     w=11 needs n>=44; sign codes A(9,3)=%d, A(10,3)=%d, A(11,3)=%d "
          "against A(8,2)=%d)" % (Abin(9, 3), Abin(10, 3), Abin(11, 3), Abin(8, 2)))
    print("      %3s %3s %10s %12s %12s %12s %s" %
          ("n", "w", "best pub.", "term", "128*A(n,8,8)", "Johnson*sgn", "conclusive?"))
    for n in range(36, 45):
        base = 128 * Acw(n, 8)
        for w in (9, 10, 11):
            if 4 * w > n:
                continue
            d = 2 * ((w + 1) // 2)
            lo = Acw(n, w)
            sgn = Abin(w, sign_dist(w))
            up = johnson_cw(n, d, w) * sgn
            print("      %3d %3d %10d %12d %12d %12d %s" %
                  (n, w, lo, lo * sgn, base, up,
                   "yes, ruled out" if up < base else "no published code is near"))
    print("    For w=10 and w=11 the Johnson bound alone is not conclusive, so here is")
    print("    the gap between what would be needed and what is published anywhere:")
    print("      %3s %3s %30s %12s %26s" %
          ("n", "w", "would need A(n,2ceil(w/2),w) >=", "Johnson <=",
           "largest published value"))
    for n in range(40, 45):
        base = 128 * Acw(n, 8)
        for w in (10, 11):
            if 4 * w > n:
                continue
            d = 2 * ((w + 1) // 2)
            sgn = Abin(w, sign_dist(w))
            need = base // sgn + 1
            doc = [(k[0], CW[k][0]) for k in CW if k[1] == d and k[2] == w
                   and "monotone" not in CW[k][1]]
            atn, biggest = max(doc, key=lambda z: z[1]) if doc else (0, 0)
            print("      %3d %3d %30d %12d %10d (published only at n=%d;"
                  " Brouwer tabulates no n>32)"
                  % (n, w, need, johnson_cw(n, d, w), biggest, atn))

    # ------------------------------------------- exact-rational geometry of the winners
    print("\n[6b] The winning chains satisfy every geometric condition exactly")
    print("     (rational arithmetic; <x,y> <= n_0/2 for all ordered level pairs):")
    from fractions import Fraction
    for n in range(32, 45):
        _, ch = results[n]
        n0 = ch[0]
        HALF = Fraction(n0, 2)
        ok = True
        for iu, wu in enumerate(ch):
            a2u = Fraction(n0, wu)                       # a_u^2 = n_0 / w_u
            # (a) same level, same support: <x,y> = a_u^2 * (w_u - 2*d_sign)
            ok &= a2u * (wu - 2 * sign_dist(wu)) <= HALF
            # (b) same level, different supports: overlap at most floor(w_u/2)
            ok &= a2u * (wu // 2) <= HALF
            # (c) different levels u < v: <x,y> <= a_u a_v w_v ; square both sides
            for iv in range(iu + 1, len(ch)):
                wv = ch[iv]
                a2v = Fraction(n0, wv)
                ok &= (a2u * a2v * wv * wv) <= HALF * HALF
        check(ok, "n=%2d chain %-12s all level-pair inner products <= n_0/2 = %s"
                  % (n, str(ch), HALF))

    # ---------------------------------------------------------------- sensitivity
    print("\n[7] Sensitivity: what would have to improve, and by how much,")
    print("    for each dimension to yield a new record.")
    print("      %3s %28s %14s %14s %10s" %
          ("n", "binding parameter", "current", "needed", "shortfall"))
    for n in range(32, 45):
        v, ch = results[n]
        need = PROJECT[n] + 1 - v
        # (a) improve A(n,8,8)
        a88 = Acw(n, 8)
        need88 = a88 + -(-need // 128) if need > 0 else a88
        # (b) improve the top code
        tops = [(Abin(n0, sign_dist(n0)), n0) for n0 in range(32, n + 1)]
        t, n0 = max(tops)
        tag = "" if need <= 1 else "   <-- ERS is not the holder here"
        print("      %3d %28s %14d %14d %10d%s" %
              (n, "A(%d,8,8)" % n, a88, need88, max(0, need), tag))
        print("      %3d %28s %14d %14d %10d" %
              (n, "A(%d,%d)  [n_0=%d]" % (n0, sign_dist(n0), n0), t,
               t + max(0, need), max(0, need)))

    print("\n[8] Consistency with the published record column:")
    for n in range(32, 45):
        v, ch = results[n]
        if n in (38, 43):
            check(v < PROJECT[n],
                  "n=%d: ERS gives %d, below the best known %d (held by a different "
                  "construction)" % (n, v, PROJECT[n]))
        else:
            check(v == PROJECT[n],
                  "n=%d: ERS with current tables gives exactly the best known %d"
                  % (n, v))

    # ------------------------------------------------------ the one structural lead
    print("")
    print("[9] The one place ERS leaves slack, and why it does not pay.")
    print("    ERS forces two level-0 supports to overlap in <= n_0/2, which is what")
    print("    makes the sign patterns of different supports unconstrained and gives the")
    print("    clean product A(n,n_0,n_0)*A(n_0,ceil(n_0/4)).  For n <= 44 and n_0 >= 30")
    print("    two n_0-subsets of [n] meet in >= 2n_0-n > n_0/2, so A(n,n_0,n_0) = 1 and")
    print("    level 0 is a SINGLE support.  One could instead take several supports with")
    print("    larger overlap and constrain the signs across them.  Worked out in the")
    print("    only case where it could matter, n = 44 with n_0 = 40:")
    print("      * two supports [44]\E_1, [44]\E_2 with E_1,E_2 disjoint 4-sets share")
    print("        36 coordinates, and <x,y> = 36 - 2*(disagreements there), so the only")
    print("        cross requirement is 'differ in >= 8 of 36' -- far weaker than the")
    print("        'differ in >= 10 of 40' required inside a single support;")
    print("      * but a 589824-word code projected to 36 coordinates saturates the cube:")
    print("        a Hamming ball of radius 7 in {0,1}^36 holds %d points, and" % vol7)
    print("        589824 * %d / 2^36 = %.0f, i.e. an average point of {0,1}^36 already"
          % (vol7, 589824.0 * vol7 / 2.0 ** 36))
    print("        has about %.0f first-support words within distance 7.  Essentially no"
          % (589824.0 * vol7 / 2.0 ** 36))
    print("        room is left for a second support of comparable size;")
    print("      * and it would take 2097152/589824 = 3.56 supports' worth merely to match")
    print("        the single-support n_0 = 44 term that is already in use.")
    print("    At n = 44 the winning chain has n_0 = n, so there is exactly one support")
    print("    and no freedom of this kind at all.")

    print("\n" + "=" * 88)
    if FAIL:
        print("FAILURES:")
        for f in FAIL:
            print("   " + f)
        return 1
    print("ALL CHECKS PASSED.")
    print("CONCLUSION: with every currently published lower bound on A(n,d,w) and")
    print("A(m,d), the Edel-Rains-Sloane construction maximised over ALL admissible")
    print("chains reproduces the best known kissing number in every dimension 32-44")
    print("and improves none of them.")
    print("=" * 88)
    return 0

if __name__ == "__main__":
    sys.exit(main())
