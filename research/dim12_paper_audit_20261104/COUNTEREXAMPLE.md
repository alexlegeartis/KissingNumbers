# A valid 48-system outside the paper's paired-floor model

This observation does **not** contradict Theorem 1 of
[Takhanov–Assylbekov–Yun, arXiv:2606.18984v2](https://arxiv.org/abs/2606.18984v2).
That theorem gives sufficient conditions. It is not a classification of all
48-systems. The distinction matters when using the paper to design searches.

Take the 48 rows of `strict48_counterexample.npy`, interpret their binary
coordinates as exact dyadic rationals, and normalize each row mathematically.
The file is copied from an archived numerical system released by the authors;
its provenance and SHA256 are in `strict48_certificate.json`.

`strict48.py` proves, entirely with integers, that the resulting unit vectors
v1,...,v48 satisfy

\[
 |v_{ik}|\leq\frac12-2^{-20},\qquad
 v_i\cdot v_j\leq\frac12-2^{-20}\quad(i\ne j).
\]

The diagnostic maximum coordinate magnitude is about 0.4999976033, and the
maximum pairwise cosine is about 0.4999964940. The rigorous bound is the
rational number 524287/1048576, independent of these floating point summaries.

For a positive dot product, the exact comparison against a rational threshold
p/q is `q^2 dot^2 <= p^2 norm_i^2 norm_j^2`. The corresponding coordinate
test is `q^2 coordinate^2 <= p^2 norm_i^2`. Different dyadic scales for different
rows cancel from these comparisons.

Theorem 1's paired floors contain, for each i, the two vectors
(1/2,b_i) and (-1/2,b_i), with squared norm of b_i equal to 3/4. Their inner
product is exactly -1/4+3/4=1/2. Every such model therefore has at least
16 contact pairs. The system certified here has **no** contact at cosine 1/2.
It cannot be isometric to any paired-floor configuration of that theorem,
including after a common rotation or a relabeling.

There is also an open set of such systems. Move each unit vector independently
by less than 2^-22 in Euclidean distance along S5. Each coordinate changes by
less than 2^-22 and each pairwise inner product by less than twice that amount.
Both inequalities stay strictly below 1/2. This is an open subset of (S5)^48,
of dimension 240. It includes independent deformations of the individual
points, rather than only the common rotations in Theorem 2.

This is a precise reason to release the paired-floor equalities and include
other archived systems in initialization. It does not by itself produce a
49-point system or a 842-point kissing arrangement.
