# An exact strict 504-point group orbit in 14D

This is a search building block, **not a new kissing-number record**: the known
14D construction contains 1,932 points. The two rational orthogonal matrices
in `discovery/recipes/group14_strict504.json` generate a 504-element subgroup
of the E7-spread group used in the repository's 1,932-point reconstruction.
The supplied integer seed has a full orbit of 504 points, exact rank 14 and
maximum normalized inner product 174,556,206 / 399,936,544 < 1/2.

From the repository root, run the independently exact orbit construction:

    python -m discovery run discovery/recipes/group14_strict504.json --out /tmp/group14-504.json
    python -m discovery verify /tmp/group14-504.json

The numerical seed search was followed by integer rationalization and exact
pair checks. The new `search.group_orbit_seeds` operator generalizes that
step to several jointly optimized free orbits. Its output is only a numerical
proposal; `groups.certify_multi_orbit` reconstructs the exact matrix group,
deduplicates the orbits and checks all relative orbit products, including
unequal seed norms. `group14_search_strict504.json` demonstrates that two-step
pipeline. The 13D and 14D multi-orbit recipes are exploratory and do not
establish any new bound.
