#!/usr/bin/env python3
"""Regenerate data/dimension27_200540.txt from the published data set (provenance).

    python scripts/rebuild_from_published.py path/to/dimensions25-31.txt [out.txt]

This script is about provenance, not about the theorem: it shows that the
configuration file was produced mechanically from Henry Cohn's published
dimension-27 block, not assembled by hand.  It reads that block, recovers the
Leech minimal vectors, the twelve cap directions, the twelve axis points and the
five 496-element classes, re-assigns four of the classes to the four triples of
one of the two admissible partitions, and writes the result.

It needs numpy and the 97 MB source file, so it is not part of the standard run;
scripts/verify_configuration.py proves the theorem without either.

Every printed token is round-tripped through the parser before being written,
and the finished file is re-read and compared against the configuration it was
built from.
"""

import os, sys, hashlib, collections, itertools as it
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import parse_token, rstr                            # noqa: E402

if len(sys.argv) < 2:
    sys.exit("usage: rebuild_from_published.py path/to/dimensions25-31.txt [out.txt]")
SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir,
    "data", "dimension27_200540.txt")

if not __debug__:
    raise SystemExit("run this script without -O: its output is checked by assertions")


def say(*a):
    print(*a, flush=True)

# --------------------------------------------------------------------------
# Z[sqrt2, sqrt3, sqrt6] in the basis (1, sqrt2, sqrt3, sqrt6)
# --------------------------------------------------------------------------
def fmul(x, y):
    a1, b1, c1, d1 = x[..., 0], x[..., 1], x[..., 2], x[..., 3]
    a2, b2, c2, d2 = y[..., 0], y[..., 1], y[..., 2], y[..., 3]
    return np.stack([a1*a2 + 2*b1*b2 + 3*c1*c2 + 6*d1*d2,
                     a1*b2 + b1*a2 + 3*(c1*d2 + d1*c2),
                     a1*c2 + c1*a2 + 2*(b1*d2 + d1*b2),
                     a1*d2 + d1*a2 + b1*c2 + c1*b2], axis=-1)

def fdot(X, Y):
    return fmul(X, Y).sum(axis=-2)

# the renderer must invert the parser on every shape it can emit
for _t in (0, 1, -1, 2, -3, 4, -6, 9, 12, -12):
    for _v in ((_t, 0, 0, 0), (0, _t, 0, 0), (0, 0, _t, 0), (0, 0, 0, _t),
               (2, _t, 0, 0), (0, 0, 2, _t), (_t, 0, 3, -1), (-2, _t, 1, 3)):
        assert parse_token(rstr(_v)) == tuple(_v), (_v, rstr(_v))
say("[render]   common.rstr round-trips against common.parse_token")

# --------------------------------------------------------------------------
# read the source block and rebuild the configuration
# --------------------------------------------------------------------------
def load_dim27(path):
    with open(path) as f:
        lines = f.readlines()
    i = next(k for k, L in enumerate(lines) if L.startswith("Dimension: 27"))
    n = int(lines[i + 1].split(":")[1])
    assert [int(x) for x in lines[i + 2].split(":")[1].strip().split(",")] == [1] * 27
    cache, R = {}, np.zeros((n, 27, 4), dtype=np.int64)
    for r, L in enumerate(lines[i + 4: i + 4 + n]):
        ps = L.strip().split(",")
        assert len(ps) == 27
        for j, p in enumerate(ps):
            if p not in cache:
                cache[p] = parse_token(p)
            R[r, j] = cache[p]
    say(f"[read]     {n} points from {path}")
    return R

R = load_dim27(SRC)
head, tail = R[:, :24, :], R[:, 24:, :]
hn, tn = fdot(head, head)[:, 0], fdot(tail, tail)[:, 0]
eq_m, cap_m, ax_m = tn == 0, (3 * tn == hn + tn) & (hn > 0), hn == 0
sc = np.zeros(len(R), dtype=np.int64)
for v, s in [(2, 12), (8, 6), (18, 4), (32, 3), (72, 2), (288, 1)]:
    sc[(hn == v) & ~ax_m] = s
for v, s in [(3, 12), (12, 6), (48, 3)]:
    sc[(tn == v) & ax_m] = s
assert (sc > 0).all()
H, T = head[..., 0] * sc[:, None], tail * sc[:, None, None]

capH, capT, AX = H[cap_m], T[cap_m], T[ax_m]
grp = collections.OrderedDict()
for i, k in enumerate(tuple(x.ravel().tolist()) for x in capT):
    grp.setdefault(k, []).append(i)
ws = np.array(list(grp), dtype=np.int64).reshape(len(grp), 3, 4)
uniq = collections.OrderedDict()
for i, v in enumerate(grp.values()):
    uniq.setdefault(frozenset(map(tuple, capH[v].tolist())), []).append(i)
S = [np.array(sorted(k), dtype=np.int64) for k in uniq]
GW = np.array([fdot(np.repeat(ws[i][None], len(ws), 0), ws) for i in range(len(ws))])[..., 0]
tri = [t for t in it.combinations(range(len(ws)), 3)
       if all(GW[a, b] == -72 for a, b in it.combinations(t, 2))]
parts = [c for c in it.combinations(tri, 4) if len(set(sum(c, ()))) == len(ws)]
assert len(parts) == 2
PART = parts[0]

Leech = np.vstack([H[eq_m]] + S)
used = set().union(*[set(map(tuple, S[i].tolist())) for i in range(4)])
EQ = np.array([r for r in Leech.tolist() if tuple(r) not in used], dtype=np.int64)
NH = np.vstack([S[g] for g, trip in enumerate(PART) for _ in trip])
NW = np.concatenate([np.repeat(ws[w][None], len(S[g]), 0)
                     for g, trip in enumerate(PART) for w in trip])
say(f"[build]    {len(AX)} axis + {len(NH)} caps + {len(EQ)} equator = "
    f"{len(AX) + len(NH) + len(EQ)}")
assert len(AX) + len(NH) + len(EQ) == 200540

# --------------------------------------------------------------------------
# rescale to the three printed shapes
# --------------------------------------------------------------------------
assert (AX % 3 == 0).all() and (EQ % 3 == 0).all()
AXo, EQo = AX // 3, EQ // 3
assert (fdot(AXo, AXo)[:, 0] == 48).all() and (fdot(AXo, AXo)[:, 1:] == 0).all()
assert ((EQo ** 2).sum(1) == 32).all()
assert ((NH ** 2).sum(1) == 288).all()
assert (fdot(NW, NW)[:, 0] == 144).all() and (fdot(NW, NW)[:, 1:] == 0).all()

rows = ([["0"] * 24 + [rstr(AXo[i, j]) for j in range(3)] for i in range(len(AXo))]
        + [[str(int(x)) for x in NH[i]] + [rstr(NW[i, j]) for j in range(3)]
           for i in range(len(NH))]
        + [[str(int(x)) for x in EQo[i]] + ["0", "0", "0"] for i in range(len(EQo))])
assert len(rows) == 200540 and all(len(r) == 27 for r in rows)

with open(OUT, "w", newline="\n") as f:
    f.write("Dimension: 27\n")
    f.write(f"Number of points: {len(rows)}\n")
    f.write("Inner product coefficients: " + ",".join(["1"] * 27) + "\n")
    f.write("Points:\n")
    for r in rows:
        f.write(",".join(r) + "\n")
toks = sorted({t for r in rows for t in r})
say(f"[write]    {OUT}: {len(rows)} points, {len(toks)} distinct coordinate tokens")

# --------------------------------------------------------------------------
# read it back and compare against what was built
# --------------------------------------------------------------------------
with open(OUT) as f:
    lines = f.read().split("\n")
assert lines[0] == "Dimension: 27"
assert lines[1] == "Number of points: 200540"
assert lines[2] == "Inner product coefficients: " + ",".join(["1"] * 27)
assert lines[3] == "Points:" and lines[-1] == "" and len(lines) == 200545
back, cache = np.zeros((200540, 27, 4), dtype=np.int64), {}
for i, L in enumerate(lines[4:200544]):
    ps = L.split(",")
    assert len(ps) == 27
    for j, p in enumerate(ps):
        if p not in cache:
            cache[p] = parse_token(p)
        back[i, j] = cache[p]
want = np.zeros((200540, 27, 4), dtype=np.int64)
want[:12, 24:] = AXo
want[12:5964, :24, 0] = NH
want[12:5964, 24:] = NW
want[5964:, :24, 0] = EQo
assert (back == want).all(), "the file does not read back as the configuration"
say("[readback] file re-parses to exactly the configuration it was built from")
say(f"[sha256]   {hashlib.sha256(open(OUT, 'rb').read()).hexdigest()}")
say(f"[tokens]   {', '.join(toks[:8])} ...")
