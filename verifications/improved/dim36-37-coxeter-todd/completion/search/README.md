# Continue the finite phase search

Python 3.9+, NumPy and a C++17 compiler suffice. Run the main package's
`verify.py` to check a certificate; these optional programs perform a heuristic
search and never publish or install a result.

Prepare a scratch directory outside the package and compile the two tools there:

```bash
pkg=/absolute/path/to/verifications/improved/dim36-37-coxeter-todd
job=/tmp/kissing-phase-refill
python "$pkg/completion/search/prepare.py" "$job" --mode equator
c++ -O3 -std=c++17 "$pkg/completion/search/filter.cpp" -o "$job/filter"
c++ -O3 -std=c++17 "$pkg/completion/search/refill.cpp" -o "$job/refill"
cd "$job"
./filter 4 6
./refill 5000000 361028 6
```

`prepare.py` defaults to the current certificate. Pass `--certificate` with the
retained `phase/certificate.json` to rebuild the original search instance. It
writes the complete pool, relative-phase dot table and incumbent in native
little-endian integer format, with the certificate hash in `pool.json`.
`filter` keeps candidates with at most four current blockers; it also writes
`expanded.bin` after any immediate insertions. `refill` writes its best selection
to `best.bin`. Its arguments specify a deterministic step count, RNG seed and
maximum allowed body dot. Steps and runtime are distinct; there is no guarantee
of improving the incumbent.

For the auxiliary latitude search use a separate directory, `--mode lifts`,
`./filter 4 2` and `./refill 5000000 361028 2`. It considers each three-point
latitude orbit separately, excludes exact conflicts with the supplied equator,
and preserves the fixed-base distance argument. A different equator requires
fresh preparation. Coupled equator/lift searches can lose latitude points, so
evaluate the **total** count when combining results.

To inspect `best.bin`, read it with
`numpy.fromfile(path,dtype=numpy.int32).reshape(-1,5)`. Install its rows in a copy
of the corresponding certificate field, update the two counts, and run both
exact checkers. The finite heuristic does not replace those checks. The compact
shipped `delta.json` is the reproducible record of the accepted trades.
