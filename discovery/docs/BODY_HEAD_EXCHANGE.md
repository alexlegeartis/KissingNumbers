# Finite body/head exchanges

`search.body_head_exchange` crosses a supplied finite body catalogue with a
supplied rational head palette. It searches direct insertions, one-for-two
exchanges and, with `max_removed: 2`, two-for-three exchanges. There is no
dimension-specific branch, target cardinality or supplied exchange in the
operator. The source bodies, head palette and family envelopes remain imported
inputs; this does not establish held-out construction discovery.

Each raw body row `x` has `x*x = body_squared_scale * body_norm_squared`.
Its geometric body is `x/sqrt(body_squared_scale)`. A head is an integer
numerator row divided by `head_denominator`. Incumbent points specify their
body indices and head numerators. The body norm is common, while the total
point norms may differ. Every input body norm and duplicate source row is
checked before searching.
The current cache backend requires an integral common body norm and integral
twice-body products against the incumbent. Unsupported rational body products
are rejected explicitly.

The optional fixed complete families have head numerators, body squared norms
and a table of upper bounds on twice the body inner product, indexed by body
class. These envelopes are supplied assumptions. Pure-head points can also
be retained. The operator does not prove family populations or that supplied
envelopes bound a geometric configuration.

Floating-point screening generates candidates and may miss feasible points.
Body/incumbent products use integer arithmetic, or exact integer sums in
binary32 under an explicit absolute-accumulation bound below 2^24. Candidate
groups share their proposed incumbent deletions. A compatible pair or triangle
can improve cardinality while charging each deleted point only once.

Before emitting an exchange, Python integers check every replacement against
all retained incumbent points, every fixed-family envelope, every pure head
and the other replacements. The sign of the dot product is retained before
squaring the normalized inequality. Arbitrarily large rational head
denominators therefore do not introduce an angular tolerance into these gates.

The result has kind `exchange_proposal` and verification status `unverified`.
Its exact local gates do not certify the supplied sources or the full output.
Independent construction replay is required before registering a lower bound.
`replacement_test_limit` is explicit; a stopped search reports truncation.
Point, product-cache, input-bit, catalogue and work budgets apply.

## Tested 43D construction

The adapter in
[`research/highdim_characters_20261104`](../../research/highdim_characters_20261104/README.md)
uses 30,240 bodies and 305 head directions from an archived incumbent. No delta
or expected count is passed to the operator. Three successive searches add
one point each. Both existing independent exact replays certify the resulting
**2,919,785-point** construction, improving the archived 2,919,782 by three.

The separate 43D task discovered the original exchange mechanism. The generic
engine independently generates the same final point set, taking the last two
exchanges in a different order. This is a seeded discovery replay. The P48n
minimum/design inputs and complete-family population identities retain the
baseline scope; large shells are not freshly enumerated.
