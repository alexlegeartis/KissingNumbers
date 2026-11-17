# Independent QR48 affine-plane augmentation audit

PASS, using the cited published QR48 minimum-distance theorem:

-88D:6,973,882,994, versus repository6,967,552,256 (+0.0908602873%).
-89D:11,622,942,176, versus9,270,596,604 (+25.3742630866%).
-90D:20,921,047,930, versus12,341,351,672 (+69.5199074301%).
-91D:20,921,062,652, versus16,435,507,572 (+27.2918561252%).

Repository values were independently read from https://github.com/alexlegeartis/KissingNumbers on2026-11-17 during the preceding audit. These are comparisons with the current displayed repository, not exhaustive worldwide-priority assertions. Code construction and simplex embedding are classical ingredients, not newly discovered parent codes.

## Parent finite checks and scope

replay_generator.py independently reruns exact factorization of x^47-1 overF3, generator rank24, self-orthogonality/self-duality, shortened ranks21 and20, and rank2 pair projection. Generated certificate hash1505ab829025e938ad09de3575656ccf5f547537d8fd5432533596fa963b22e8 exactly matches the producer. verify_sparse.py additionally checks irreducibility of the degree23 factor. The computed Gleason weight enumerator is conditional on the known minimum theorem, not an independent minimum computation.

Primary source: G.Nebe, On Extremal Self-Dual Ternary Codes of Length48, International Journal of Combinatorics2012,154281, https://doi.org/10.1155/2012/154281 and https://arxiv.org/html/1109.1681 . Introduction and section3 establish actual known extremal ternary self-dual[48,24,15] codes; Lemma3.5 identifies the extended QR realization used here. Minimum15 is inherited from this published theorem.

## Dense layers

Shortening3 or4 zero coordinates gives[45,21,>=15] and[44,20,>=15]. Exact dimensions also follow because any dependency on those coordinates would be a dual word of weight<=4, contrary to minimum15. Unit equilateral-triangle embedding of each symbol has norm squared N=45 or44 and scalar product N-3h/2. The range15<=h<=N gives absolute scalar product<=N/2. Global negatives therefore double the sets without duplication: no triangle vertex is the negative of another triangle vertex.

For89D, in the45-block construction impose equality of the first Cartesian coordinate in the first two blocks. Five of nine symbol pairs satisfy it. Uniformity follows from full projection rank onto the three shortened and two selected positions, or the verified shortened pair rank. Hence the antipodal count is10*3^19. Rotate the two equal coordinates into one coordinate (x1+x2)/sqrt2 and discard their zero difference. All retained coordinates have absolute value<=sqrt2. For88,90,91 this same bound is valid, although<=1 suffices.91D embeds the90D dense layer with an extra zero coordinate.

Every dense vector has at least44 nonzero coordinates: each original triangle block contributes at least one;89D merging may remove one, leaving at least44. Thus no dense vector equals any sparse vector below, whose support has at most4.

## Sparse affine-plane layer

Label n Cartesian coordinate axes by distinct integers0,...,n-1 viewed asF2^7 vectors. Select every four-element support with label XOR zero. Each three-element subset determines its unique fourth label; therefore distinct supports intersect in at most2 axes. Put all16 sign choices on each support, each nonzero coordinate having magnitude sqrt(N)/2, and also add the2n signed axis vectors of magnitude sqrt(N).

Every vector has squared norm N. Two sign patterns on the same support differ in at least one sign and have inner product<=N/2. Distinct supports contribute at most2*(N/4)=N/2. Axis/support products are at mostN/2; axis/axis products are0 or-N. All16 signs already include each antipode exactly once. Different supports cannot duplicate, and axis points are disjoint from four-support points.

A dense-sparse cross product is at most4*(sqrt(N)/2)*sqrt2=2 sqrt(2N)<=N/2 because N>=32. Dense-axis cross products are at mostsqrt(2N)<=N/2. These bounds apply to both signs. Therefore the union is a valid kissing configuration with no duplicates.

## Independent counting

verify_sparse.py uses a different count from producer triple completion. It groups coordinate pairs by their XOR. Distinct pairs of equal XOR are disjoint; every valid four-support appears in exactly three pair partitions. Thus support count is(1/3)*sum_x binomial(number_of_pairs_with_XOR_x,2). A second direct enumeration of all four-tuples and XOR test agrees. Counts are19,626;20,458;21,334;22,254 supports in88–91D. Sparse totals16*supports+2n are314,192;327,506;341,524;356,246. All16-sign same-support inequalities were checked exactly, as were rational squared cross bounds and total arithmetic. sparse_report.json records results.
