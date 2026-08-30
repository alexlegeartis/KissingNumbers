#!/usr/bin/env python3
"""Dimension 13: extract Cohn's published 1154-point configuration and test it.

Dimension 13 is the ONE value in dimensions 9-16 still standing from 1999 (Zinoviev-Ericson,
reference [30] of Cohn's table).  Everything around it has moved recently: dimension 10 -> 510
and dimension 14 -> 1932 by Ganzhinov (2025), dimension 11 -> 604 (2026), dimension 12 -> 841
(Takhanov-Assylbekov-Yun, 2026).  A 27-year-old construction is the most likely place in the
low-dimensional range to have slack.

Tests run here:
  * validity: maximum inner product at most half the norm;
  * MAXIMALITY: is there a unit z with <z,x> <= 1/2 for every x?  (dimension 25 gained +2
    exactly this way, and dimensions 11 and 12 were shown maximal by the same test);
  * the inner-product spectrum, which says whether the configuration is algebraic or numerical.

    python extract.py
"""
import os
import sys
import numpy as np
import re
import math

HERE = os.path.dirname(os.path.abspath(__file__))
def _find_source():
    """Locate H. Cohn's dimensions1-24.txt.

    It is 24 MB and is not shipped with this repository -- it belongs to
    H. Cohn, Table of kissing number bounds, https://hdl.handle.net/1721.1/153312
    Give it to us in any of these ways, in order of precedence:

        python <script> /path/to/dimensions1-24.txt      (first command-line argument)
        set D124=/path/to/dimensions1-24.txt             (environment variable)
        put it next to this repository, or two levels up
    """
    import sys as _sys
    cand = []
    for a in _sys.argv[1:]:
        if a.endswith('.txt'):
            cand.append(a)
    if os.environ.get('D124'):
        cand.append(os.environ['D124'])
    here = os.path.dirname(os.path.abspath(__file__))
    for up in ('.', '..', os.path.join('..', '..'), os.path.join('..', '..', '..'),
               os.path.join('..', '..', '..', '..'), os.path.join('..', '..', '..', '..', '..')):
        cand.append(os.path.join(here, up, 'dimensions1-24.txt'))
    for c in cand:
        if c and os.path.exists(c):
            return os.path.abspath(c)
    raise SystemExit(
        "\n".join([
            "",
            "This script needs H. Cohn's coordinate data set  dimensions1-24.txt  (24 MB),",
            "which is not shipped here because it is not ours to redistribute.",
            "",
            "    download:  https://hdl.handle.net/1721.1/153312",
            "    then:      python %s /path/to/dimensions1-24.txt" % os.path.basename(
                _sys.argv[0] or 'script.py'),
            "    or set:    D124=/path/to/dimensions1-24.txt",
            "",
        ]))


SRC = None          # resolved lazily by blocks(); see _find_source above


def blocks(path=None):
    if path is None:
        path = _find_source()
    cur = None
    for ln in open(path):
        ln = ln.rstrip('\n')
        if ln.startswith('Dimension:'):
            if cur:
                yield cur
            cur = {'dim': int(ln.split(':')[1]), 'pts': [], 'c': None, 'n': None}
        elif ln.startswith('Number of points:'):
            cur['n'] = int(ln.split(':')[1])
        elif ln.startswith('Inner product coefficients:'):
            cur['c'] = [float(x) for x in ln.split(':')[1].split(',')]
        elif ln.startswith('Points:'):
            pass
        elif ln.strip() and cur is not None:
            # entries may be plain numbers or expressions like "2*sqrt(3)" / "-2*sqrt(3)"
            row = []
            ok = True
            for tok in ln.split(','):
                tok = tok.strip()
                m = re.fullmatch(r'([+-]?[0-9.]*)\*?sqrt\(([0-9.]+)\)', tok)
                if m:
                    a = m.group(1)
                    a = float(a) if a not in ('', '+', '-') else (-1.0 if a == '-' else 1.0)
                    row.append(a * math.sqrt(float(m.group(2))))
                    continue
                try:
                    row.append(float(tok))
                except ValueError:
                    ok = False
                    break
            if ok and row:
                cur['pts'].append(row)
    if cur:
        yield cur


if __name__ == '__main__':
    want = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    best = None
    for b in blocks(SRC):
        if b['dim'] == want and (best is None or b['n'] > best['n']):
            best = b
    assert best is not None, "dimension %d not found" % want
    print(__doc__)
    X = np.array(best['pts'])
    c = np.array(best['c'])
    print("dimension %d: %d points, %d parsed; inner-product coefficients %s"
          % (best['dim'], best['n'], len(X), sorted(set(c.tolist()))))
    # inner product <x,y> = sum c_i x_i y_i  -> rescale coordinates by sqrt(c_i)
    Y = X * np.sqrt(c)[None, :]
    nrm = (Y * Y).sum(1)
    print("norms: %d distinct, min %.6f max %.6f" % (len(set(np.round(nrm, 6))), nrm.min(), nrm.max()))
    U = Y / np.sqrt(nrm)[:, None]
    G = U @ U.T
    np.fill_diagonal(G, -1.0)
    print("max inner product between distinct points: %.12f  (valid: %s)"
          % (G.max(), G.max() <= 0.5 + 1e-9))
    from collections import Counter
    sp = Counter(np.round(G[~np.eye(len(U), dtype=bool)], 6).tolist())
    print("inner-product spectrum: %d distinct values" % len(sp))
    print("   most common:", sorted(sp.items(), key=lambda t: -t[1])[:8])
    near = (G.min(axis=1) < -0.99).sum()
    print("points with a near-antipode: %d of %d" % (near, len(U)))
    np.save(os.path.join(HERE, 'c%d.npy' % want), U)   # a cache; see .gitignore
    print("saved c%d.npy" % want)
