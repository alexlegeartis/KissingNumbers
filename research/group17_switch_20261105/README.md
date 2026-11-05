# An inequivalent 5,730-point configuration in 17D

This is a **neutral 128-for-128 exchange**, not a new lower bound. Both the
source and switched configuration have 5,730 points, matching the known 17D
lower bound. The switched set differs from the source in 128 points, and the
two configurations are not isometric: their overall pairwise inner-product
histograms agree, but their distributions of per-point inner-product profiles
differ. This certificate establishes inequivalence to the included source;
it does not claim inequivalence to every configuration in the literature.

`coordinates.npz` contains integer source and switched rows for the metric
`diag(1,...,1,2)`; every row has squared norm 93,312 and every distinct
inner product is at most 46,656. Divide transformed coordinates by the square
root of 93,312 after identifying the weighted metric with Euclidean space.

The switched points were selected individually from the source and seven
rational isometric images. In the weighted metric, use reflections with
normals e_0+e_16, e_1+e_16, and e_2+e_3+e_4+e_5. The finite words used are
[], [0], [1], [2], [0,1], [0,2], [1,2], [0,1,2], applied left to right. The
normals mix the weighted coordinate into the other coordinates. Whole orbits
were not required: the search charged shared deletions and chose individual
points from the image union.

Run `python verify.py` in this directory. It independently checks every norm
and every pair for both codes in exact integer arithmetic, checks all 128
replacements belong to the stated group-image union, and verifies the
per-point-profile inequivalence certificate. `generate_pool.py` rebuilds the
finite image graph, and `python search.py 17w 170705` runs an individual-point
shared-deletion search. The search is heuristic and a failed run gives no upper
bound. The related `search.rational_image_exchanges` operator is registered in
the main discovery engine.
