#!/usr/bin/env python3
"""The best known lower bound on tau(n) for every n from 48 to 96, with its source.

Cohn's table has 48 and 72 and NOTHING in between or after.  This fills the gap from the
three places the answer actually lives: this repository's 47 claims, the direct-sum /
monotonicity floors of common/published.py where nothing is claimed, and -- at n = 96 --
Edel-Rains-Sloane, which takes over there.

    python table48-96.py            # markdown table on stdout
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DELIV = os.path.join(HERE, '..', '..', '..')
sys.path.insert(0, os.path.join(DELIV, 'common'))
import published as PUB

# ---- this repository's claims, read from the generated RESULTS.md (never typed by hand)
CLAIM, ROW = {}, {}
for line in open(os.path.join(DELIV, 'RESULTS.md'), encoding='utf-8'):
    m = re.match(r'\| (\d+) \| ([\d\s\u202f\u00a0]+) \| \*\*([\d\s\u202f\u00a0]+)\*\* \|', line)
    if m:
        d = int(m.group(1))
        CLAIM[d] = int(re.sub(r'\D', '', m.group(3)))
        ROW[d] = re.search(r'\[([^\]]+)\]\(', line).group(1) if '[' in line else ''

# ---- the one dimension above the claims
ERS96 = 12886999232
# computed by cap96.py, which prints it; kept here as a constant because that
# script is a driver with nothing to import.  It went stale once, when the class
# family grew and dimension 95 rose past it, so check it against cap96.py after
# any change to the class families.
CAP96 = 6480558568

print("| n | tau(n) >= | source |")
print("|---|---|---|")
run = 0
for n in range(48, 97):
    if n == 96:
        v, src = ERS96, ("**Edel-Rains-Sloane**, chain (96,24,6,1), with codetables.de "
                         "2026-08-26 -- [96,33,24] gives A(96,24) >= 2^33 on its own")
    elif n in CLAIM:
        v, src = CLAIM[n], "**this work** -- %s" % ROW[n]
    else:
        v, src = PUB.floor_for(n)
    assert v >= run, "monotonicity fails at n = %d" % n
    run = v
    print("| %d | %s | %s |" % (n, "{:,}".format(v).replace(",", " "), src))

print()
print("dimension 96 detail")
print("  Gamma_72 (+) Leech, the direct-sum floor        6 218 372 160")
print("  the cap construction over Gamma_72 at k = 24    %s   (this repository's"
      % "{:,}".format(CAP96).replace(",", " "))
print("                                                                  mechanism, and NOT")
print("                                                                  a record)")
print("  ERS level 0 alone, the [96,33,24] code          8 589 934 592")
print("  ERS chain (96,24,6,1)                          %s   <- the entry" % "12 886 999 232")
