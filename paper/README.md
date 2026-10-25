# `paper/` — the write-up, and the three checkers that hold it to the repository

[`kissing46.tex`](kissing46.tex) — *New lower bounds for kissing numbers in dimensions 18
through 96*. 52 pages, `amsart`, no packages beyond `amssymb`, `amsmath`, `booktabs`,
`longtable` and `hyperref`. References are [`kissing46.bib`](kissing46.bib), formatted with
`amsplain`, so a full build is

```bash
pdflatex kissing46 && bibtex kissing46 && pdflatex kissing46 && pdflatex kissing46
```

[`kissing46.bbl`](kissing46.bbl) is tracked as well, although it is generated: it is what a
build without `bibtex` reads, and without it one `pdflatex` pass prints `[?]` for all
forty-five citations. [`kissing46.pdf`](kissing46.pdf) is a build of the `.tex` beside it, kept for convenience.
Nothing here checks a PDF against its source — rebuild with the command above if you want
to confirm it; the build is 52 pages, with no undefined citation or reference.
The generated current-results table and the historical appendix produce a few
typographic over/underfull-box warnings; they do not affect the mathematical
content.

## Current results and historical appendix

The main paper now reports the 63 registered entries that remain state of the
art in the reviewed sources. `current_table.tex` is generated from the current
65-row `RESULTS.md` registry and includes the October advances in dimensions
26--30, 38--39 and 64--68, along with the 28D, 36--37D and 41--46D updates. The 25 October
73D construction combines the completed half-coset with a verified non-lattice
cap replacement and a norm-ten addition, reaching 6,218,204,076. Exact
whole-shell norm-ten existence proofs also improve dimensions 77--79. The two superseded current
entries remain in the paper's historical appendix with their construction details
and the frozen September table. The new midpoint certificates reach
216,568 in 29D and 224,280 in 30D; the latitude/root replacements reach
605,346 in 38D and 768,924 in 39D. Their proofs and verification scope
are included in the current construction appendix.

`values.tex` continues to pin the September numerical snapshot because the
historical derivations and their checks depend on it. `current_values.tex` and
`current_table.tex` are generated from `RESULTS.md` by `mkcurrentvalues.py` and
`mkcurrenttable.py`; run both scripts with `--check` to detect drift.

## One bound, one place: [`values.tex`](values.tex)

Every claimed bound used to be typed into `kissing46.tex` four times — the introduction, the
theorem, its row of Table 4, and the count identity that derives it — with its gain in two more
tables. Moving one dimension meant six edits in five places, and the only thing between that and
a wrong number in print was `factcheck.py` noticing afterwards. Now the manuscript writes

```latex
$K(31)\ge\Kh{31}$          % 238 662, this work
$\Kp{31}$                   % 238 350, previously published
$\Kgain{31}$                % +312
```

and [`values.tex`](values.tex) defines all of them, **generated** by
[`mkvalues.py`](mkvalues.py) from [`registry_snapshot.json`](registry_snapshot.json)
for the manuscript’s 53 dimensions, and from [`../RESULTS.md`](../RESULTS.md) for
new dimensions. Updating the manuscript also requires updating its snapshot
and derivations. 93 sites in the
manuscript now come from that one file (32 `K(d)≥` statements, 52 rows of Table 4, the net-gain
row of Table 2 and the largest-gain cell of every row of Table 1), and `python mkvalues.py
--check` exits non-zero if it has drifted from those inputs.

`factcheck.py` **expands** these macros before it parses, so all of its checks still see the
numbers a reader sees — macroising the manuscript cannot quietly turn a check into a no-op — and
it fails if any `\Kh{...}` is left unexpanded or if `values.tex` disagrees with the frozen snapshot and current registry.
The addends of a count identity (`175\,904+61\,968+118+672`) stay literal: those are a
derivation, not a bound.

## The checkers

Three of these are jobs of [`../run_all.py`](../run_all.py) and exit non-zero on a mismatch;
the fourth reports counts and never fails. Run them from this directory.

```bash
python factcheck.py      # every figure in the paper, recomputed from the repository
python unsupported.py    # no figure in the paper is without a home in the repository
python formulas.py       # every formula in the paper, re-derived with sympy
python style_scan.py -v  # the prose faults listed in the writing guidelines
```

`factcheck.py` (774 checks) does not parse the repository's conclusions, it **re-derives**
them: the moment identities, the count of the layered construction, Caro–Wei, the
Edel–Rains–Sloane totals, the bound on independent sets, and every entry of Tables 1, 3 and 4.
`formulas.py` (257 checks) checks the *algebra* rather than the numbers — the moment identity
on all four shells, the
seven inequalities of Proposition 5.1, the threshold ⌊μ/3⌋/μ, the direction bound
γ/(2γ−1), weak duality on random instances, and the Ozeki affine relation. **Run all four
after any edit to the paper.**

Several guards are relational rather than literal, because a guard that pins a number in
place makes the number un-correctable: where a sentence states two quantities and a
consequence, the consequence is recomputed from the sentence's own figures.

## One term differs between the paper and the code

The paper calls a set of minimal lines pairwise at |cos| ≤ γ an **independent set**, in the
graph Γ_γ(L) it defines in §5.1. The repository calls the same object a **class**, in
directory names, script names and variables. They are the same thing.

## Where this sits

The paper cites the repository it is in — the ordinary arrangement for a paper with a
verification package. It lived in a sibling directory of the author's working tree until
2026-08-30, which meant these checkers skipped in every clone; moving it here is what lets a
reader run them. `../audit.py` reads `kissing46.tex` too, for two of its documentation
checks.

`registry_snapshot.json` pins the historical manuscript's 53 numerical results.
The current registry is authoritative for the main table. `mkvalues.py`,
`mkcurrentvalues.py`, `mkcurrenttable.py`, `factcheck.py` and `unsupported.py`
check the two layers and require the current registry to retain or improve every
historical result.
