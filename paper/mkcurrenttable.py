#!/usr/bin/env python3
"""Generate the main current-SOTA results table from RESULTS.md."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mkcurrentvalues import read_rows, tex_number  # noqa: E402

OUT = os.path.join(HERE, "current_table.tex")

# All 66 current claims exceed their reviewed comparison. Earlier mechanisms
# remain in the historical appendix, including dimensions regained this month.
CONSTRUCTIONS = {
    40: r"three affine lattice shells (\S\ref{ss:d40current})",
    31: r"midpoints and dual E7 centres (\S\ref{ss:d31current})",
    25: r"collective rational deformation (\S\ref{ss:d25current})",
    18: r"bent-coset hexagon over odd $\mathrm{BW}_{16}$ (\S\ref{sec:d18})",
    26: r"shared midpoint replacements (\S\ref{ss:d2627current})",
    27: r"shared midpoint replacements (\S\ref{ss:d2627current})",
    28: r"joint Leech layers (\S\ref{ss:d28current})",
    29: r"shared midpoint replacements (\S\ref{ss:d2930current})",
    30: r"shared midpoint replacements (\S\ref{ss:d2930current})",
    38: r"joint latitude families (\S\ref{ss:d3839current})",
    39: r"deformed sign orbits (\S\ref{ss:d3839current})",
    36: r"Coxeter--Todd phase code (\S\ref{ss:d3637})",
    37: r"Coxeter--Todd phase code with latitude layers (\S\ref{ss:d3637})",
    41: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
    42: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
    43: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
    44: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
    45: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
    46: r"$P_{48n}$ heads (\S\ref{ss:p48current})",
}
CONSTRUCTIONS.update({d: r"layers over $P_{48}$ (\S\ref{ss:big})"
                      for d in range(49, 61)})
CONSTRUCTIONS.update({d: r"binary code layers (\S\ref{ss:codecurrent})"
                      for d in (61, 62, 63, 66, 67, 68, 95, 96)})
CONSTRUCTIONS.update({d: r"signed-support code layer (\S\ref{ss:codecurrent})"
                      for d in (64, 65)})
CONSTRUCTIONS.update({d: r"$\Gamma_{72}$ cross-section (\S\ref{ss:g72})"
                      for d in (69, 70, 71)})
CONSTRUCTIONS.update({d: r"nearest $\Gamma_{72}$ cosets (\S\ref{ss:cosetscurrent})"
                      for d in range(73, 87)})
CONSTRUCTIONS.update({d: r"$\Gamma_{72}$ cap layers (\S\ref{ss:big})"
                      for d in range(87, 95)})

CONSTRUCTIONS[73] = r"$\Gamma_{72}$ norm-ten addition (\S\ref{ss:cosetscurrent})"

CONSTRUCTIONS.update({d: r"$\Gamma_{72}$ norm-ten existence bound (\S\ref{ss:cosetscurrent})"
                      for d in (77, 78, 79)})

SOTA = sorted(CONSTRUCTIONS)


def row(dimension, previous, current):
    return (r"$%d$ & $\Kcurprev{%d}$ & $\Kcur{%d}$ & %s \\" %
            (dimension, dimension, dimension, CONSTRUCTIONS[dimension]))


def render():
    rows = read_rows()
    missing = sorted(set(SOTA) - set(rows))
    if missing:
        raise SystemExit("current registry is missing: %s" % missing)
    lines = [
        r"\scriptsize",
        r"\begin{longtable}{rrrp{0.30\textwidth}}",
        r"\caption{Current SOTA lower bounds in the reviewed sources.  The",
        r"reviewed set contains 66 certified claims; earlier constructions are retained in",
        r"Appendix~\ref{app:historical}.}\label{tab:all}\\",
        r"\toprule",
        r"$n$ & recorded comparison & this work & construction\\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"$n$ & recorded comparison & this work & construction\\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for d in SOTA:
        previous, current = rows[d]
        lines.append(row(d, previous, current))
    lines.append(r"\end{longtable}")
    lines.append(r"\normalsize")
    return "\n".join(lines) + "\n"


def main():
    text = render()
    if "--check" in sys.argv:
        with open(OUT, encoding="utf-8") as fh:
            current = fh.read()
        if current == text:
            print("current_table.tex is up to date")
            return 0
        print("current_table.tex is stale; run: python paper/mkcurrenttable.py")
        return 1
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("wrote current_table.tex: %d rows" % len(SOTA))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
