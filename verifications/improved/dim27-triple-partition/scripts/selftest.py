#!/usr/bin/env python3
"""Negative controls: check that the verifier rejects configurations that are wrong.

    python scripts/selftest.py

A checker that always says PASS proves nothing.  This script takes the real
200540-point file, confirms it is accepted, and then makes one targeted
corruption at a time, confirming each is rejected.  Every corruption replaces a
row rather than adding one, so the file keeps its declared size and the failure
is the specific defect being tested rather than a size mismatch.

The corruptions are the ways this family of configurations can actually break:

  * a repeated point, and the subtler case of the same point written at a
    different scale, which a textual duplicate check would miss;
  * condition (E) broken -- a Leech vector used as a cap head is left on the
    equator, putting an equatorial point at inner product sqrt(2/3) from a cap;
  * a cap moved to a direction outside its triple, so two caps carrying one
    class sit at cosine +1/2 instead of -1/2;
  * a foreign vector inserted into a cap class, breaking the cosine 1/4 bound
    that caps on a common direction must satisfy;
  * a head of the right shape and norm whose support is not an octad, so it is
    not a Leech minimal vector;
  * every axis point moved onto a cap direction, which is condition (A) at its
    worst: sqrt3 lies in the coordinate ring, so t/sqrt3 has squared norm 48 and
    integral coordinates, and the resulting configuration has 5952 pairs at
    inner product 1/sqrt3 while every other check still passes.  An adversarial
    review found that this was the one case the checker got wrong, so it is now
    the control that matters most;
  * a header whose declared count disagrees with the file;
  * a coordinate whose spelling is ambiguous -- an unsigned juxtaposition such
    as "2sqrt(3)", which could be read as a sum or as a product, and a token
    written with a non-ASCII digit.  Both were accepted by an earlier version of
    the parser, which would have let the file mean one thing to this package and
    another to any other reader.

Each case runs the checker on a full-size file, so this takes a couple of
minutes.  No third-party packages are required.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import parse_token, rstr                                  # noqa: E402
ROOT = os.path.join(HERE, os.pardir)
DATA = os.path.join(ROOT, "data", "dimension27_200540.txt")
CHECKER = os.path.join(HERE, "verify_configuration.py")
TMP = os.path.join(ROOT, "selftest-tmp")
N_AXIS, N_CAPS = 12, 5952


if not __debug__:
    raise SystemExit("run this script without -O: its fixtures are checked by assertions")


def write(name, header, rows):
    os.makedirs(TMP, exist_ok=True)
    p = os.path.join(TMP, name)
    with open(p, "w", newline="\n", encoding="utf-8") as f:
        f.write("\n".join(header) + "\n")
        f.write("\n".join(rows) + "\n")
    return p


def run(path):
    r = subprocess.run([sys.executable, CHECKER, path], capture_output=True, text=True)
    out = [l.strip() for l in r.stdout.split("\n") if l.strip()]
    verdict = ""
    for l in out:
        if "checks passed" in l or l.startswith("VERIFIED"):
            verdict = l
    if not verdict:
        tail = (out or [""])[-1]
        err = [l.strip() for l in r.stderr.split("\n") if l.strip()]
        verdict = (err[-1] if err else tail)
    return r.returncode, verdict[:96]


def main():
    if not os.path.exists(DATA):
        sys.exit(f"missing {DATA}")
    text = open(DATA).read().rstrip("\n").split("\n")
    header, rows = text[:4], text[4:]
    if len(rows) != 200540:
        sys.exit(f"unexpected data file: {len(rows)} rows")
    axis = list(range(N_AXIS))
    caps = list(range(N_AXIS, N_AXIS + N_CAPS))
    equator = list(range(N_AXIS + N_CAPS, len(rows)))

    def integer_row(i):
        return all(t.lstrip("+-").isdigit() for t in rows[i].split(","))

    assert integer_row(equator[0]), "equator rows should be all-integer"

    cases = [("the file as published", header, rows, True)]

    # 1. a repeated point
    d = list(rows)
    d[equator[5]] = d[equator[9]]
    cases.append(("a repeated point", header, d, False))

    # 2. the same point at a different scale (rows are normalised, so still a repeat)
    d = list(rows)
    d[equator[5]] = ",".join(str(3 * int(t)) for t in rows[equator[9]].split(","))
    cases.append(("the same point written at three times the scale", header, d, False))

    # 3. condition (E): leave a cap head on the equator
    caphead = rows[caps[0]].split(",")[:24]
    assert all(int(x) % 3 == 0 for x in caphead)
    d = list(rows)
    d[equator[7]] = ",".join([str(int(x) // 3) for x in caphead] + ["0", "0", "0"])
    cases.append(("condition (E) broken: a cap head left on the equator", header, d, False))

    # 4. a cap moved to a direction outside its triple
    c0 = rows[caps[0]].split(",")
    foreign = next(rows[i].split(",")[24:] for i in caps
                   if rows[i].split(",")[24:] != c0[24:])
    d = list(rows)
    d[caps[1]] = ",".join(rows[caps[1]].split(",")[:24] + foreign)
    cases.append(("a cap moved to a foreign direction", header, d, False))

    # 5. a foreign vector inserted into a cap class
    victim = rows[equator[0]].split(",")[:24]
    d = list(rows)
    d[caps[2]] = ",".join([str(3 * int(x)) for x in victim] + c0[24:])
    cases.append(("a foreign vector inserted into a cap class", header, d, False))

    # 6. a head of shape (+-2^8) whose support is not an octad.  The checker
    #    derives its Golay code as the span of the supports of the shape-(+-2^8)
    #    heads, so "not an octad" means "not one of the supports occurring in
    #    the file".  Assert that rather than trusting the name of the case: were
    #    positions 0..7 an octad, this control would be testing nothing.
    supports = set()
    for line in rows:
        cs = line.split(",")[:24]
        if not all(t.lstrip("+-").isdigit() for t in cs):
            continue
        vals = [abs(int(t)) for t in cs]
        nz = [i for i, v in enumerate(vals) if v]
        if len(nz) == 8 and len({vals[i] for i in nz}) == 1:
            supports.add(frozenset(nz))
    assert len(supports) == 759, f"expected 759 octads in the file, found {len(supports)}"
    assert frozenset(range(8)) not in supports, "positions 0..7 are an octad after all"
    non_octad = ["2"] * 8 + ["0"] * 16
    d = list(rows)
    d[equator[11]] = ",".join(non_octad + ["0", "0", "0"])
    cases.append(("an equator head whose support is not an octad", header, d, False))

    # 7. every axis point moved onto a cap direction.  A cap tail t has squared
    #    norm 144; t/sqrt3 has squared norm 48, and since sqrt3 is in the ring
    #    it has integral coordinates whenever 3 divides the rational and sqrt2
    #    parts of every coordinate of t.  Then t.(t/sqrt3) = 144/sqrt3 = 48 sqrt3,
    #    far above the true bound of 72, so those caps sit at cosine 1/sqrt3.
    def over_sqrt3(tok):
        a, b, c, d = parse_token(tok)          # a + b r2 + c r3 + d r6
        if a % 3 or b % 3:                     # times sqrt3/3
            return None
        return rstr((c, d, a // 3, b // 3))

    tails = []
    seen = set()
    for i in caps:
        t = tuple(rows[i].split(",")[24:])
        if t not in seen:
            seen.add(t)
            tails.append(t)
    moved = [tuple(over_sqrt3(x) for x in t) for t in tails]
    if len(tails) == 12 and all(all(x is not None for x in m) for m in moved):
        d = list(rows)
        for k in range(12):
            d[axis[k]] = ",".join(["0"] * 24 + list(moved[k]))
        cases.append(("every axis point moved onto a cap direction", header, d, False))
    else:
        print("  [skip] axis-on-cap-direction case could not be constructed")

    # 8. header count disagrees with the file
    bad_header = list(header)
    bad_header[1] = "Number of points: 200541"
    cases.append(("header count disagrees with the file", bad_header, rows, False))

    # 9-10. the coordinate syntax itself.  "2sqrt(3)" has two readings, 2 + sqrt3
    #    and 2*sqrt3, and a token spelled with a non-ASCII digit is a number to a
    #    Unicode-aware \d but not to a human or to another tool.  Either would
    #    let this file mean one thing here and something else elsewhere, so the
    #    parser has to refuse both rather than silently pick a reading.
    d = list(rows)
    d[equator[13]] = ",".join(["2sqrt(3)"] + rows[equator[13]].split(",")[1:])
    cases.append(("a coordinate written as an unsigned juxtaposition, '2sqrt(3)'",
                  header, d, False))
    d = list(rows)
    d[equator[15]] = ",".join(["\uff12"] + rows[equator[15]].split(",")[1:])
    cases.append(("a coordinate written with a non-ASCII digit, U+FF12",
                  header, d, False))

    print("=" * 78)
    print("negative controls for scripts/verify_configuration.py")
    print("=" * 78, flush=True)
    wrong = 0
    for i, (name, hdr, rs, should_pass) in enumerate(cases):
        p = write(f"case{i:02d}.txt", hdr, rs)
        code, verdict = run(p)
        os.remove(p)
        accepted = (code == 0)
        ok = accepted == should_pass
        wrong += 0 if ok else 1
        print(f"  [{'ok  ' if ok else 'WRONG'}] {name}")
        print(f"         expected {'accept' if should_pass else 'reject'}, "
              f"got {'accept' if accepted else 'reject'}", flush=True)
        print(f"         {verdict}", flush=True)
    wrong += bridge_controls(header, rows)

    print()
    if wrong == 0:
        print("All controls behaved correctly: the checkers accept the published")
        print("configuration and reject every one of these corruptions.")
    else:
        print(f"{wrong} CONTROL(S) WRONG.")
    try:
        os.rmdir(TMP)
    except OSError:
        pass
    return 1 if wrong else 0


def bridge_controls(header, rows):
    """Negative controls for verify_paper_claims.py.

    These are the exploits an adversarial review produced against an earlier
    version of that script, which checked the equator only by cardinality and
    rescaled rows without guards.  Each must now be rejected.
    """
    import json
    checker = os.path.join(HERE, "verify_paper_claims.py")
    jsrc = os.path.join(ROOT, "data", "construction.json")
    os.makedirs(TMP, exist_ok=True)
    doc = json.load(open(jsrc))

    def run2(conf, js):
        r = subprocess.run([sys.executable, checker, conf, js],
                           capture_output=True, text=True)
        out = [l.strip() for l in r.stdout.split("\n") if l.strip()]
        v = ""
        for l in out:
            if "checks passed" in l:
                v = l
        if not v:
            err = [l.strip() for l in r.stderr.split("\n") if l.strip()]
            v = err[-1] if err else (out or [""])[-1]
        return r.returncode, v[:96]

    def data(name, rs, hdr=None):
        p = os.path.join(TMP, name)
        with open(p, "w", newline="\n", encoding="utf-8") as f:
            f.write("\n".join(hdr or header) + "\n")
            f.write("\n".join(rs) + "\n")
        return p

    cases = []
    # the axis points of construction.json forged onto the cap directions
    forged, ok = [], True
    for d in doc["cap_directions"]:
        coords = []
        for a, b, c, dd in d:
            if a % 3 or b % 3:
                ok = False
            coords.append([c, dd, a // 3, b // 3])
        forged.append(coords)
    if ok:
        alt = dict(doc)
        alt["axis_points"] = forged
        jp = os.path.join(TMP, "forged.json")
        json.dump(alt, open(jp, "w", newline="\n"), separators=(",", ":"))
        cases.append(("construction.json axis points forged onto cap directions",
                      data("b_clean.txt", rows), jp))
    # a vector of squared norm 32 that is not in the Leech lattice, on the equator
    d = list(rows)
    d[N_AXIS + N_CAPS + 3] = ",".join(["5", "2", "1", "1", "1"] + ["0"] * 19 + ["0"] * 3)
    cases.append(("a non-Leech vector of norm 32 on the equator",
                  data("b_leech.txt", d), jsrc))
    # an equator row whose squared norm is not a valid rescaling
    d = list(rows)
    d[N_AXIS + N_CAPS + 5] = ",".join(["4", "4", "1"] + ["0"] * 21 + ["0"] * 3)
    cases.append(("an equator row of squared norm 33",
                  data("b_norm.txt", d), jsrc))
    # duplicated rows with an inflated header
    hdr = list(header)
    hdr[1] = "Number of points: 200545"
    cases.append(("5 duplicated rows, header says 200545",
                  data("b_dup.txt", rows + [rows[-1]] * 5, hdr), jsrc))

    print()
    print("negative controls for scripts/verify_paper_claims.py")
    print("-" * 78, flush=True)
    wrong = 0
    for name, conf, js in cases:
        code, verdict = run2(conf, js)
        ok = code != 0
        wrong += 0 if ok else 1
        print(f"  [{'ok  ' if ok else 'WRONG'}] {name}")
        print(f"         expected reject, got {'reject' if ok else 'accept'}")
        print(f"         {verdict}", flush=True)
    for f in os.listdir(TMP):
        os.remove(os.path.join(TMP, f))
    return wrong


if __name__ == "__main__":
    sys.exit(main())
