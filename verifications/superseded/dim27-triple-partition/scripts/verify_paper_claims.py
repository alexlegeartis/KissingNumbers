#!/usr/bin/env python3
"""Check the finite assertions of the paper, and bridge them to the data file.

    python scripts/verify_paper_claims.py [config.txt [construction.json]]

This script is secondary to scripts/verify_configuration.py, which proves the
theorem on its own.  What this one adds is the link between the paper's small
combinatorial data -- the twelve cap directions, the twelve axis points, the
five 496-element classes and the way classes are assigned to directions -- and
the 200540 rows of data/dimension27_200540.txt.

It reads data/construction.json and checks, in exact arithmetic and with no
third-party packages:

  * the twelve cap directions form a cuboctahedron, and so do the twelve axis
    points, and the two share no vertex;
  * exactly eight triples of directions are pairwise at cosine -1/2, no four
    directions are pairwise at cosine at most -1/2, and those eight triples
    admit exactly two partitions of the twelve directions;
  * the five classes are pairwise disjoint of size 496, with internal cosines
    at most 1/4 and cross cosines reaching exactly 1/2 (the hypotheses of the
    optimality theorem);
  * the published assignment groups the directions as {3,3,2,2,2}, giving
    m = 7 and 196560 + 7*496 + 12 = 200044;
  * the new assignment groups them as {3,3,3,3}, is one of the two partitions,
    gives m = 8 and 196560 + 8*496 + 12 = 200540;
  * among cap-cap and axis-axis pairs, those attaining inner product exactly
    1/2 are the three cap families of 900096, 277408 and 11904 ordered pairs and
    48 ordered axis pairs (equator-equator pairs attain 1/2 as well, in far
    greater numbers, and are not counted here);
  * and finally that rebuilding the configuration from this data reproduces
    data/dimension27_200540.txt exactly.
"""
import collections
import itertools as it
import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (Report, parse_configuration, rdot, rmul, is_rational,
                    sha256sum, split_families, certify_leech_minimal)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, os.pardir, "data")
JSON = os.path.join(DATA, "construction.json")
CONF = os.path.join(DATA, "dimension27_200540.txt")


def ring(v):
    return tuple(int(x) for x in v)


def gram(vs):
    return [[rdot(a, b) for b in vs] for a in vs]


def profile(G, scale):
    """How many rows of a rational Gram see each cosine profile.

    One profile shared by every row is what vertex-transitivity looks like from
    the Gram alone; it exhibits no automorphism, and is reported as such.  The
    vectors here are triples, so they span at most 3 dimensions by construction
    and no rank check is needed.
    """
    out = collections.Counter()
    for row in G:
        c = collections.Counter(Fraction(x[0], scale) for x in row)
        out[tuple(sorted(c.items()))] += 1
    return out


def histogram(A, B, same):
    """Counts of head dot products between two lists of Leech vectors."""
    h = collections.Counter()
    if same:
        for i in range(len(A)):
            u = A[i]
            for j in range(i + 1, len(A)):
                h[sum(map(int.__mul__, u, A[j]))] += 1
    else:
        for u in A:
            for v in B:
                h[sum(map(int.__mul__, u, v))] += 1
    return h


def main(conf=None, jsonp=None):
    conf = conf or CONF
    jsonp = jsonp or JSON
    t0 = time.time()
    rep = Report("the paper's finite claims, and the bridge to the data file")

    rep.section("input")
    rep.note(f"{os.path.basename(jsonp)}  sha256 {sha256sum(jsonp)}")
    doc = json.load(open(jsonp))
    dirs = [tuple(ring(c) for c in d) for d in doc["cap_directions"]]
    axis = [tuple(ring(c) for c in a) for a in doc["axis_points"]]
    classes = [[tuple(int(x) for x in v) for v in c] for c in doc["classes"]]
    pub_assign = doc["published_assignment"]
    new_assign = doc["new_assignment"]
    kept = doc["classes_kept"]
    dropped = doc["class_dropped"]

    # -------------------------------------------------------- basic shapes
    rep.section("the small data")
    rep.check(len(dirs) == 12 and all(rdot(d, d) == (144, 0, 0, 0) for d in dirs),
              "12 cap directions, each of squared norm 144")
    rep.check(len(axis) == 12 and all(rdot(a, a) == (48, 0, 0, 0) for a in axis),
              "12 axis points, each of squared norm 48")
    rep.check(len(classes) == 5 and all(len(c) == 496 for c in classes),
              "5 classes of 496 vectors")
    rep.check(all(sum(x * x for x in v) == 32 for c in classes for v in c),
              "every class vector has squared norm 32")
    sets = [set(c) for c in classes]
    rep.check(all(sets[i].isdisjoint(sets[j]) for i in range(5) for j in range(i + 1, 5)),
              "the five classes are pairwise disjoint")
    rep.check(len(set().union(*sets)) == 5 * 496,
              "their union has 2480 vectors")

    # ------------------------------------------------------- cuboctahedra
    rep.section("the two cuboctahedra")
    GD = gram(dirs)
    rep.check(all(is_rational(x) for row in GD for x in row),
              "the cap-direction Gram matrix is rational")
    CUBO = {((Fraction(-1), 1), (Fraction(-1, 2), 4), (Fraction(0), 2),
             (Fraction(1, 2), 4), (Fraction(1), 1))}
    pD = profile(GD, 144)
    rep.check(set(pD) == CUBO and list(pD.values()) == [12],
              "all 12 cap directions see the same cuboctahedron cosine profile",
              "-1 once, -1/2 four times, 0 twice, 1/2 four times, 1 once")
    GA = gram(axis)
    rep.check(all(is_rational(x) for row in GA for x in row),
              "the axis Gram matrix is rational")
    pA = profile(GA, 48)
    rep.check(set(pA) == CUBO and list(pA.values()) == [12],
              "so do all 12 axis points: a second cuboctahedron")
    # A cap direction has squared norm 144 and an axis point 48, so d = k a
    # forces k^2 = 3, that is k = sqrt3 -- which lies in this ring.  The
    # coincidence is therefore possible and has to be tested, not argued away.
    SQRT3 = (0, 0, 1, 0)
    MSQRT3 = (0, 0, -1, 0)
    same = sum(1 for d in dirs for a in axis for k in (SQRT3, MSQRT3)
               if all(dc == rmul(ac, k) for dc, ac in zip(d, a)))
    rep.check(same == 0, "no axis point lies on a cap direction, either sign",
              "a = t/sqrt3 has squared norm 48 and integral coordinates, so this is real")

    # ------------------------------------------------------- triple search
    rep.section("triples and partitions")
    tri = [t for t in it.combinations(range(12), 3)
           if all(GD[a][b] == (-72, 0, 0, 0) for a, b in it.combinations(t, 2))]
    quad = [q for q in it.combinations(range(12), 4)
            if all(GD[a][b][0] <= -72 for a, b in it.combinations(q, 2))]
    parts = [c for c in it.combinations(tri, 4) if len(set(sum(c, ()))) == 12]
    rep.check(len(tri) == 8, "exactly 8 triples are pairwise at cosine -1/2")
    rep.check(len(quad) == 0, "no 4 directions are pairwise at cosine at most -1/2")
    rep.check(len(parts) == 2, "those 8 triples admit exactly 2 partitions of the 12")
    rep.check(sorted(map(list, tri)) == sorted(doc["triples_at_cos_minus_half"]),
              "construction.json lists the same 8 triples")
    rep.check(sorted(sorted(map(list, p)) for p in parts)
              == sorted(sorted(map(list, p)) for p in doc["partitions_into_four_triples"]),
              "and the same 2 partitions")

    # ---------------------------------------------------- the two groupings
    rep.section("the published grouping and the new one")
    pub_groups = collections.defaultdict(list)
    for d, c in enumerate(pub_assign):
        pub_groups[c].append(d)
    rep.check(sorted((len(v) for v in pub_groups.values()), reverse=True) == [3, 3, 2, 2, 2],
              "the published assignment groups the 12 directions as {3,3,2,2,2}")
    m_pub = 12 - len(pub_groups)
    rep.check(m_pub == 7, "so m = 12 - 5 = 7", str(m_pub))
    rep.check(196560 + 7 * 496 + 12 == 200044,
              "and 196560 + 7*496 + 12 = 200044", str(196560 + 7 * 496 + 12))

    new_groups = collections.defaultdict(list)
    for d, c in enumerate(new_assign):
        new_groups[c].append(d)
    rep.check(sorted((len(v) for v in new_groups.values()), reverse=True) == [3, 3, 3, 3],
              "the new assignment groups them as {3,3,3,3}")
    rep.check(sorted(tuple(sorted(v)) for v in new_groups.values())
              in [sorted(tuple(sorted(t)) for t in p) for p in parts],
              "and that grouping is one of the two partitions")
    m_new = 12 - len(new_groups)
    rep.check(m_new == 8, "so m = 12 - 4 = 8", str(m_new))
    rep.check(196560 + 8 * 496 + 12 == 200540,
              "and 196560 + 8*496 + 12 = 200540", str(196560 + 8 * 496 + 12))
    rep.check(len(kept) == 4 and dropped not in kept and len(set(kept)) == 4,
              "four distinct classes are kept and the fifth is dropped",
              f"kept {kept}, dropped {dropped}")
    rep.check((196560 + 8 * 496 + 12) - (196560 + 7 * 496 + 12) == 496,
              "the gain is exactly 496, the deletions the dropped class was costing")
    rep.check(all(GD[a][b][0] <= -72
                  for ds in pub_groups.values() for a, b in it.combinations(ds, 2)),
              "in the published assignment, directions sharing a class are at cosine <= -1/2")
    rep.check(all(GD[a][b] == (-72, 0, 0, 0)
                  for ds in new_groups.values() for a, b in it.combinations(ds, 2)),
              "in the new assignment they are at cosine exactly -1/2")

    # --------------------------------------------- class Grams, hypotheses
    rep.section("the hypotheses of the optimality theorem")
    hist_self = [histogram(classes[i], None, True) for i in range(5)]
    hist_cross = {}
    for i in range(5):
        for j in range(i + 1, 5):
            h = histogram(classes[i], classes[j], False)
            hist_cross[(i, j)] = hist_cross[(j, i)] = h
    rep.check(all(max(h) == 8 for h in hist_self),
              "inside every class the maximum cosine is exactly 1/4",
              f"maxima {[max(h) for h in hist_self]} over 32")
    rep.check(all(max(h) == 16 for h in hist_cross.values()),
              "every pair of distinct classes attains cosine exactly 1/2",
              "this is the cross-pair hypothesis of Theorem 4(ii)")

    # --------------------------------------------- pairs attaining 1/2
    rep.section("cap-cap and axis-axis pairs attaining inner product exactly 1/2")
    # new_assign indexes the four kept classes; map into the published numbering
    cls_of_dir = {d: kept[new_assign[d]] for d in range(12)}
    same_dir = sum(2 * hist_self[cls_of_dir[d]][8] for d in range(12))
    half_dir = 0
    minus_dir = 0
    shared_at_plus_half = 0
    for a in range(12):
        for b in range(12):
            if a == b:
                continue
            ca, cb = cls_of_dir[a], cls_of_dir[b]
            g = GD[a][b][0]
            if g == 72:                       # cosine 1/2
                if ca == cb:
                    # two directions at cosine +1/2 carrying one class would put
                    # a cap at 2/3 + 1/6 = 5/6 against itself on the other
                    # direction; Lemma 2(i) forbids it.  Count it rather than
                    # passing over it, so the enumeration below cannot be
                    # silently incomplete.
                    shared_at_plus_half += 1
                else:
                    half_dir += hist_cross[(ca, cb)][16]
            elif g == -72:                    # cosine -1/2
                if ca == cb:
                    minus_dir += 496
    rep.check(shared_at_plus_half == 0,
              "no class is carried by two directions at cosine +1/2 (Lemma 2(i))",
              f"{shared_at_plus_half} such ordered pairs")
    rep.check(same_dir == 900096,
              "<w,w'> = 1 and <u,u'> = 1/4: 900096 ordered pairs", str(same_dir))
    rep.check(half_dir == 277408,
              "<w,w'> = 1/2 and <u,u'> = 1/2: 277408 ordered pairs", str(half_dir))
    rep.check(minus_dir == 11904,
              "<w,w'> = -1/2 and u = u': 11904 ordered pairs",
              f"{minus_dir} = 4 classes x 6 ordered pairs x 496")
    n_axis = sum(1 for i in range(12) for j in range(12)
                 if i != j and GA[i][j] == (24, 0, 0, 0))
    rep.check(n_axis == 48, "48 ordered axis pairs at cosine 1/2", str(n_axis))

    # ----------------------------------------------------------- the bridge
    rep.section("bridge: rebuild the configuration from this data")
    cfg = parse_configuration(conf, expect_dimension=27)
    # the same guarded splitter verify_configuration.py uses: exact rescaling,
    # integrality, and the 1/3 ratio that identifies a cap row
    try:
        eq_list, cap_list, ax_list, _, _ = split_families(cfg)
    except ValueError as e:
        rep.check(False, "the file splits into equator, cap and axis families", str(e))
        return rep.done("")
    rep.check(len(cfg.rows) == 200540, "the file has 200540 rows", str(len(cfg.rows)))
    rep.check((len(eq_list), len(cap_list), len(ax_list)) == (194576, 5952, 12),
              "split as 194576 equator + 5952 caps + 12 axis",
              f"{len(eq_list)} + {len(cap_list)} + {len(ax_list)}")
    file_eq, file_caps, file_ax = set(eq_list), set(cap_list), set(map(tuple, ax_list))
    rep.check(len(file_eq) == len(eq_list) and len(file_caps) == len(cap_list)
              and len(file_ax) == len(ax_list),
              "no row is repeated within its family")

    want_caps = set()
    for c, ds in new_groups.items():
        for d in ds:
            for y in classes[kept[c]]:
                want_caps.add((y, dirs[d]))
    rep.check(file_caps == want_caps,
              "the 5952 caps are exactly the prescribed (class, direction) pairs",
              f"{len(file_caps)} in the file, {len(want_caps)} prescribed")
    rep.check(file_ax == set(map(tuple, axis)),
              "the 12 axis points are exactly the prescribed ones")
    used = set().union(*[set(classes[kept[c]]) for c in new_groups])
    rep.check(len(used) == 1984, "the four used classes contribute 1984 deleted vectors")
    rep.check(file_eq.isdisjoint(used),
              "no used class vector remains on the equator")
    rep.check(len(file_eq) + len(used) == 196560,
              "equator and deleted together number 196560",
              f"{len(file_eq)} + {len(used)}")
    rep.check(set(classes[dropped]) <= file_eq,
              "the dropped class is back on the equator: its 496 vectors are equatorial")
    # the count above is only a count; certify that those 196560 really are the
    # minimal vectors of a copy of Lambda_24, so that "the whole set" is earned
    cert = certify_leech_minimal(sorted(file_eq | used))
    rep.check(cert["norms_32"] and cert["distinct"],
              "all 196560 are distinct integer vectors of squared norm 32")
    rep.check(cert["census_ok"], "their shape census is 1104 + 97152 + 98304")
    rep.check(cert["n_octads"] == 759 and cert["golay_dim"] == 12
              and cert["weights_ok"] and cert["octads_are_weight8"],
              "their octad supports form the extended binary Golay code")
    rep.check(cert["all_in_leech"],
              "and every one satisfies the congruences defining Lambda_24",
              "so equator + deleted really is a full set of Leech minimal vectors")

    print()
    print(f"  elapsed {time.time() - t0:.1f}s")
    return rep.done(f"The paper's finite claims hold, and {os.path.basename(conf)} is\n"
                    "           exactly the configuration they describe.")


if __name__ == "__main__":
    a = sys.argv[1:]
    sys.exit(main(a[0] if a else None, a[1] if len(a) > 1 else None))
