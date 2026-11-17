# Ternary outer-code words with nonlinear hexagonal phase layers

This package proves ordinary kissing-number lower bounds **21,943,640,934 in90D**, **21,943,640,936 in91D**, and **21,943,640,940 in92D**. It replaces the preceding sparse layers; their cross compatibility is not assumed and their points are not added. No global optimality or exhaustive worldwide-priority claim is made.

| Dimension | Attributable external-source baseline | New bound | External source |
|---:|---:|---:|---|
| 90 | 6,218,183,254 | **21,943,640,934** | [Nebe Gamma72](https://arxiv.org/abs/1008.2862v3) plus [18D tail](../../../common/published.py), derived orthogonal sum |
| 91 | 6,218,187,548 | **21,943,640,936** | [Nebe Gamma72](https://arxiv.org/abs/1008.2862v3) plus [19D tail](../../../common/published.py), derived orthogonal sum |
| 92 | 8,589,934,592 | **21,943,640,940** | [Duursma–Greferath–Litsyn–Schmidt nonlinear Golay code](https://doi.org/10.1109/18.923742), shortened |

These comparisons use the repository's attributable pre-existing external-source baselines, not our own preceding constructions. A derived sum or shortening is identified as such, not represented as an independently published record in that exact dimension.

## Construction and source premises

The dense90D layer has2·3^21=20,920,706,406 points: shorten three coordinates of the classical ternary QR48 code, embed the remaining45 symbols as unit equilateral triangles, and include global negatives. Its published minimum15 premise and independently verified generator are retained in the [preceding package](../dim68-91-lattice-code-continuation/).

Take130,944 distinct weight15 words from the ternary Pless code shortened in coordinates0,1,2. Membership, weight, distinctness and shortening zeros are explicitly checked; shell completeness is unnecessary. The systematic generator is exactly the Paley conference construction [I|S] overF3, with S indexed by infinity andF23, first row1, first column−1, and S_ij the Legendre symbol of i−j. The Pless[48,24,15] minimum is a published premise: [Nebe, On extremal self-dual ternary codes of length48](https://doi.org/10.1155/2012/154281), [author manuscript](https://arxiv.org/html/1109.1681). Rank and self-duality alone are not used as a proof of minimum15.

At each occupied block use one of three roots t0−t2,t1−t0,t2−t1, with sign determined by the outer ternary word. Each root has squared norm3. Assign phases from the explicit nonlinear ternary(15,7812,5) code of [Laaksonen and Östergård](https://www.cs.helsinki.fi/u/ahslaaks/codes/), [original data archive](https://www.cs.helsinki.fi/u/ahslaaks/codes/codes2a.zip), member3_15_5_7812.txt. Both checkers verify all30,509,766 phase pairs exactly. This supplies130,944·7,812=1,022,934,528 sparse points.

Different outer words have squared distance at least3·15=45; different phase words over the same outer word have squared distance at least9·5=45. Every dense/sparse block inner product has absolute value at most3/2, so15 occupied blocks give at most45/2. All points have norm squared45. Zero sparse blocks distinguish them from dense words. In91D add two orthogonal poles; in92D add a regular hexagon in the orthogonal two-dimensional complement, contributing six points, not six antipodal pairs.

## Immutable evidence and replay

Original production and independent packages are preserved unchanged. Their research-stage notes saying “pending” or “unsubmitted” are historical; the later independent nonlinear-code PASS is in `hexagon_sparse15_independent_20261117/AUDIT.md`. Original audit notes may compare intermediate internal constructions, but the comparison table above deliberately uses external sources.

Run `python verify.py` (Python3.10+, NumPy and SymPy). It invokes unchanged `hexagon_sparse15_independent_20261117/verify.py`, `hexagon_sparse15_phase7812_20261117/verify.py`, and `hexagon_sparse15_independent_20261117/verify_phase7812.py`, checks the exact Pless generator, replays the retained QR48 generator, and verifies the orthogonal complement and totals. The optional production regeneration script is `hexagon_sparse15_20261117/build.py`; all input code data and original checksums are shipped. Further body lengths included in the original data are finite supporting artifacts, not additional registered claims in this release.

`final.py` supplies the three declared totals to the root registry; it is not a replacement for the proof replay.

This is a retained construction; the direct latitude package supersedes its91D and92D totals.
