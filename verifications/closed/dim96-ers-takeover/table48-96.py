#!/usr/bin/env python3
"""Registered lower bounds and recorded floors in dimensions 48 through 96.

The October live Cohn page stops at 48. This table combines the current
repository registry with the recorded direct-sum and monotonicity floors of
common/published.py. It retains the old ERS and cap calculations below the table.

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

# ---- retained historical construction values
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
    if n in CLAIM:
        v, src = CLAIM[n], "**this work** -- %s" % ROW[n]
    else:
        v, src = PUB.floor_for(n)
    assert v >= run, "monotonicity fails at n = %d" % n
    run = v
    print("| %d | %s | %s |" % (n, "{:,}".format(v).replace(",", " "), src))

print()
print("## Retained dimension-96 comparisons")
print()
print("```text")
print("  Gamma_72 (+) Leech, the direct-sum floor        6 218 372 160")
print("  the cap construction over Gamma_72 at k = 24    %s   (this repository's"
      % "{:,}".format(CAP96).replace(",", " "))
print("                                                                  mechanism, and NOT")
print("                                                                  a record)")
print("  ERS level 0 alone, the [96,33,24] code          8 589 934 592")
print("  ERS chain (96,24,6,1), retained history        %s" % "12 886 999 232")
print("```")
