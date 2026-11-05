# Spectral representation discovery

`search.spectral_graph` treats the vertices of a supplied symmetric graph as a permutation representation. It computes integer or single-real-quadratic-field eigenspaces and their exact spectral projectors, searches unions of a bounded number of sectors, groups vertices with equal projection norm, and optionally adds antipodes. A finite independent-set search removes pairs exceeding the kissing threshold. Every returned Gram matrix is independently checked for positive semidefiniteness, exact rank, common norm and all angular inequalities.

The included test constructs the triangular graph on the 28 two-subsets of an 8-set. Its seven-dimensional eigenspace yields 56 antipodal points with maximum inner product 1/3. This is a sanity check, not a new lower bound.

A second sanity check constructs the regular pentagon from the graph C5 in the real quadratic field Q(sqrt(5)). Its projector certificate verifies exact symmetry, idempotence, trace and every pairwise inequality. Its reported dimension is the ambient projector rank; the selected points may span a smaller subspace. This is not a new lower bound.

This backend currently requires the graph spectrum to be integral or contained in one real quadratic field. It neither automatically discovers the graph or group nor handles general algebraic-number eigenspaces, and matrix size is limited by the configured work budgets. The numerical optimizer is only a proposal mechanism: successful output is checked exactly, and failure does not prove optimality.

The tests also provide the icosahedral graph and recover the known 12-point configuration in dimension 3 over Q(sqrt(5)). `e8_benchmark.py` builds the known E8 contact graph from its known roots and recovers a certified 240-point code in dimension 8. These are supplied-graph reconstruction benchmarks, not blind discoveries and not new records. Integer eigenvalue proposals may be estimated numerically for speed, but their annihilating polynomial and projector identities are checked in exact integer arithmetic before use.
