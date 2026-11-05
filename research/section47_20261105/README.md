# 47D sections and projected layers from P48

This is a restricted calculation, **not a new kissing-number bound**. Let S be
the norm-six minimal shell of an extremal even unimodular 48-dimensional
lattice. The inputs |S| = 52,416,000, integral inner products, minimum six,
and its spherical 11-design property are standard lattice facts; see
[Nossek, Spherical designs and lattices](https://arxiv.org/abs/1111.0772) and
[Boyvalenkov–Cherkashin](https://arxiv.org/abs/2312.05121). The latter gives
the known 47D lower bound 23,766,960 via a minimal-vector hyperplane.

Choose minimal vectors u,v with u·v = -2 and let w = u+v, of norm eight.
The verification reconstructs the joint distribution of (x·u,x·v) using
all 11-design moments, excluding ±u and ±v. It has one free parameter tau,
but two nonnegative multiplicities are tau and -6 tau, so tau=0. Hence
exactly 20,582,240 shell vectors lie in w-perp, and each x·w = ±1 layer has
12,816,384 vectors. Orthogonal projection of each layer to w-perp, followed
by its own normalization, is a kissing configuration.

Only original inner product 3 creates a conflict across the central and ±1
projected layers. `duals.json` gives rational dual polynomials. `verify.py`
checks every coefficient inequality exactly and derives degree intervals for
these bipartite conflict graphs. A matching argument bounds the union of the
+1 and -1 layers alone by 15,601,840, below the existing 47D record.

For all three layers together, the same degree inequalities only yield an
upper bound of 24,286,647 for this **restricted projected-shell pool**. This
is not a construction and does not prove that the pool can reach it. It
leaves 519,687 points of theoretical headroom above the known 47D lower bound;
actual joint selection and more precise interaction structure remain open.
Other layers, different normals, other lattices and all general configurations
are outside this bound. The full-type Cauchy search was insufficient to beat
the central layer and is not claimed exhaustive beyond that sufficient test.

Run `python research/section47_20261105/verify.py`. It needs SymPy and checks
all moment and dual identities with rational arithmetic; no numerical LP
status is used for the proof.
