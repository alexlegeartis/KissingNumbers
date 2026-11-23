# Grouped four-site designs in 65–66 dimensions

Replace the four-site sign layer of the [65D atlas](../dim65-68-complete-hyperplane-atlas/) and the [66D affine layer](../dim66-68-nested-affine-face-trades/). All other vectors stay fixed.

| Dimension | Old bound | Four-site supports | Added points | New bound |
|---|---:|---:|---:|---:|
| 65 | 418,995,074 | 10,726 | 4,960 | **419,000,034** |
| 66 | 429,284,740 | 11,418 | 5,280 | **429,290,020** |

Make eleven groups of six points. Develop the thirty literal starter blocks on eleven groups of two points cyclically. Inflate each point threefold: the four new labels lie in Z3 and sum to zero, so every triple using three different groups occurs once. Within each pair of six-point groups, take all pairs of edges from corresponding factors of a five-factor one-factorization of K6. This covers every triple using exactly two groups once. Put three more quadruples within each group. They leave precisely eight triples per group. Delete the last coordinate for 65D.

Use all sixteen signs on each support. Distinct supports meet in at most two places; two distinct points on the same support have raw dot at most two. Against every unchanged old family the worst possible dot follows from its entries: at most four for squared norms 16 and 64, two for a coordinate pole, and 5 sqrt(2) for a hyperplane head of squared norm 72.

The 66D four-site layer is optimum in its alphabet. At each point the graph of uncovered triples through it has even degrees and edge count 1 mod 3, hence at least four edges. Count across all 66 points to obtain at least 88 uncovered triples, exactly the leave in this construction. A global random sign assignment retains one sixteenth of any signed code in expectation and leaves an unsigned packing, extending the cap to arbitrary partial sign selections. This is not an upper bound on the unrestricted kissing number.

Run `python3 verify.py` to replay the retained layers and the new design. Use `--construction-only` to check the new design and arithmetic without replaying the two predecessors.
