# Independent88D three-symbol D4 compression audit

PASS: tau88 >=9,298,405,928, comprising9,298,091,736 dense points and314,192 sparse points. Frozen independent input SHA2563439267ca1dbb7d43149760db292042c47a9e4cbe4f1dbf7f050d88290c1a4d5. The twelve-pattern witness is saved as input.json, not read from a potentially mutable search output during replay.

## Source and population

Parent is the published ternary extended quadratic-residue[48,24,15] code, cited in G.Nebe2012, DOI10.1155/2012/154281, https://arxiv.org/html/1109.1681 . The preceding independent audit reproduced generator hash1505ab829025e938ad09de3575656ccf5f547537d8fd5432533596fa963b22e8, factorization, irreducibility, self-duality and ranks. Minimum15 remains a published premise, not a new exhaustive computation.

The first three original coordinates are shortened to zero. The next three coordinates specify the head prefix. Independent modular elimination confirms original generator rank24 and projection onto all six positions rank6. Therefore each prescribed prefix after shortening has exactly3^(24-6)=3^18=387,420,489 codewords. Twelve distinct prefix patterns give12*3^18 words before antipodal completion. No assumption of an uncomputed weight enumerator is used.

## Geometry

The remaining42 ternary coordinates are embedded as unit triangles, giving an84-dimensional body of squared norm42. Each selected prefix is assigned the recorded D4 root r of squared norm2, scaled by sqrt(3/2), giving a four-dimensional head of squared norm3. The resulting88-dimensional vectors all have squared norm45.

For distinct prefix patterns at Hamming distance d, original code distance implies body Hamming distance h>=15-d. The body scalar product is42-3h/2. The witness satisfies -1<=r.s<=2-d, independently checked for all66 prefix pairs. Thus total scalar product is at most42-(3/2)(15-d)+(3/2)(2-d)=45/2, and at least-21-3/2=-45/2. The verifier checks every integer h between15-d and42. For two distinct words in the same prefix fiber, h>=15 and the common head contributes3, yielding total45-3h/2 within[-45/2,45/2].

Global antipodal completion is therefore valid and gives exactly2*12*3^18=9,298,091,736 points. No duplicate or antipodal collision occurs: the42-body triangle alphabet has no antipodal alphabet pair, so a mapped positive word cannot equal the negative of another. Within a prefix, distinct bodies remain distinct; different prefixes have distinct recorded roots. The pair-angle bounds also rule out any coincident distinct source points.

## Sparse augmentation

Use exactly the preceding audited88-axis XOR-zero four-support layer:19,626 supports,16 signs per support, plus176 axis points, now scaled to norm squared45. Independent XOR-pair partition counting again gives314,192 points. Within-layer angles remain valid under scaling. Each dense body coordinate has absolute value<=1 and each head coordinate<=sqrt(3/2)<=sqrt2. Therefore dense/four-support scalar products are at most2sqrt(90)<=45/2, and dense/axis products at mostsqrt(90)<=45/2. Every dense vector has at least42 nonzero body coordinates, whereas sparse supports have at most4; hence the layers are disjoint. All sparse antipodes are already included exactly once by their sign labels.

## Comparison and scope

The independently read live repository README on2026-11-17 lists88D6,967,552,256. This construction exceeds it by2,330,853,672 (about33.453%). It also supersedes the immediately preceding audited ternary affine-plane88D result6,973,882,994. No exhaustive worldwide-priority conclusion follows from those comparisons. The new ingredient under audit is the finite twelve-prefix D4 head witness, not a new parent code.

Files verify.py and report.json provide exact finite rank, pair geometry, population and augmentation checks. The support-layer theorem and primary code source are also documented in the preceding ternary_qr48_independent_20261117/AUDIT.md.
