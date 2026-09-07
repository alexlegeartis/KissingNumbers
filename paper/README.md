# `paper/` — the write-up, and the three checkers that hold it to the repository

[`kissing46.tex`](kissing46.tex) — *New lower bounds for kissing numbers in dimensions 25
through 96*. 37 pages, `amsart`, no packages beyond `amssymb`, `amsmath`, `booktabs`,
`longtable` and `hyperref`. References are [`kissing46.bib`](kissing46.bib), formatted with
`amsplain`, so a full build is

```bash
pdflatex kissing46 && bibtex kissing46 && pdflatex kissing46 && pdflatex kissing46
```

[`kissing46.bbl`](kissing46.bbl) is tracked as well, although it is generated: it is what a
build without `bibtex` reads, and without it one `pdflatex` pass prints `[?]` for all
forty-five citations. [`kissing46.pdf`](kissing46.pdf) is a build of the `.tex` beside it, kept for convenience.
Nothing here checks a PDF against its source — rebuild with the command above if you want
to confirm it; the build is 37 pages with no warnings.

## The checkers

Three of these are jobs of [`../run_all.py`](../run_all.py) and exit non-zero on a mismatch;
the fourth reports counts and never fails. Run them from this directory.

```bash
python factcheck.py      # every figure in the paper, recomputed from the repository
python unsupported.py    # no figure in the paper is without a home in the repository
python formulas.py       # every formula in the paper, re-derived with sympy
python style_scan.py -v  # the prose faults listed in the writing guidelines
```

`factcheck.py` (665 checks) does not parse the repository's conclusions, it **re-derives**
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
