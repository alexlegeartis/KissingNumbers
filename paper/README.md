# `paper/` — the write-up, and the three checkers that hold it to the repository

[`kissing46.tex`](kissing46.tex) — *New lower bounds for kissing numbers in dimensions 18
through 96*. 47 pages, `amsart`, no packages beyond `amssymb`, `amsmath`, `booktabs`,
`longtable` and `hyperref`. References are [`kissing46.bib`](kissing46.bib), formatted with
`amsplain`, so a full build is

```bash
pdflatex kissing46 && bibtex kissing46 && pdflatex kissing46 && pdflatex kissing46
```

[`kissing46.bbl`](kissing46.bbl) is tracked as well, although it is generated: it is what a
build without `bibtex` reads, and without it one `pdflatex` pass prints `[?]` for all
forty-five citations. [`kissing46.pdf`](kissing46.pdf) is a build of the `.tex` beside it, kept for convenience.
Nothing here checks a PDF against its source — rebuild with the command above if you want
to confirm it; the build is 47 pages, with no LaTeX warnings and no undefined
citation or reference (9 over/underfull boxes, which is typography, not content).

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
[`mkvalues.py`](mkvalues.py) from [`../RESULTS.md`](../RESULTS.md), which `audit.py` in turn
generates from the packages. So updating a bound is: re-run `audit.py`, run `python
mkvalues.py`, rebuild — no prose is touched unless the *mechanism* changed. 93 sites in the
manuscript now come from that one file (32 `K(d)≥` statements, 52 rows of Table 4, the net-gain
row of Table 2 and the largest-gain cell of every row of Table 1), and `python mkvalues.py
--check` exits non-zero if it has drifted from `RESULTS.md`.

`factcheck.py` **expands** these macros before it parses, so all of its checks still see the
numbers a reader sees — macroising the manuscript cannot quietly turn a check into a no-op — and
it fails if any `\Kh{...}` is left unexpanded or if `values.tex` disagrees with `RESULTS.md`.
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

`factcheck.py` (770 checks) does not parse the repository's conclusions, it **re-derives**
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
