# Independent D3/D5 ternary head-compression audit

PASS: tau87>=6,973,870,208 and tau89>=12,397,783,154.

Frozen87 input SHA25656bc28f5641628611a0246b22f574c70ecd15dbe78c1763eb1273ec9e8c5720c;89 input SHA2567f4b94d67c5896ae9b39c4b6ebc6cfcc5aa12551547956bfa2f8d5019638c1b6. Immutable independent copies input87.json and input89.json match. verify.py and report.json preserve the exact replay.

## Parent and fibers

Use the previously independently checked generator of the published extremal ternary QR48 code, with inherited minimum distance15. Source: G.Nebe2012, DOI10.1155/2012/154281, https://arxiv.org/html/1109.1681 . The generator factorization/rank/self-duality audit remains distinct from the cited minimum theorem; there is no new exhaustive minimum computation.

Independent modular elimination verifies generator rank24 and projection onto its first six coordinates rank6. In87D shorten the first3 coordinates and prescribe the next3; in89D shorten the first2 and prescribe the next4. In either case each selected prefix fiber contains exactly3^18=387,420,489 words. There are9 and16 distinct selected prefixes respectively.

Self-orthogonality implies every codeword has weight divisible by3 because each nonzero ternary entry squares to1. Thus the distance between any two original codewords is divisible by3; shortening zero positions does not change this. This divisibility is essential for the improved lower scalar-product condition in the present witnesses.

## Head/body angle proof

Both constructions retain42 triangle body blocks, squared norm42 and dimension84. Replace the prefix by sqrt(3/2) times its recorded D3 or D5 root. The head has squared norm3 and dimension3 or5, so total dimension87 or89 and squared norm45.

For two distinct source words with prefix distance d and body distance h, h+d>=15, h<=42, and h+d is divisible by3. Body scalar product is42-3h/2; head scalar product is(3/2) times integer root dot product z. The sufficient exact conditions are z<=2-d and z>=hmax-43, where hmax=3 floor((42+d)/3)-d. The verifier checks these and every allowed integer h separately, including same-prefix distinct-word cases. All45 unordered pair types for87D and136 for89D pass, giving total inner product between-45/2 and45/2.

For a word and its own negative the scalar product is-45. Adding all global negatives therefore preserves the kissing bound. No antipodal collisions occur because the42-body triangle alphabet contains no antipodal alphabet pair. Distinct words in one prefix have distinct bodies; distinct prefixes have distinct roots. Dense counts are2*9*3^18=6,973,568,802 and2*16*3^18=12,397,455,648.

## Sparse layer

All head coordinates have absolute value<=sqrt(3/2), body coordinates<=1. The preceding audited XOR-zero four-support layer and signed axes therefore apply with the conservative coordinate bound sqrt2 and norm squared45. Dense/sparse scalar product is at most2sqrt90<45/2; dense/axis bound is smaller. Each dense vector has at least42 nonzero body coordinates, so it cannot coincide with a support-size<=4 sparse vector. All16 signs on each four-support already include each antipode once; axes are separately disjoint.

Independent XOR-pair partition counting and direct four-tuple XOR enumeration agree:18,827 supports in87D and20,458 in89D. Sparse counts16*supports+2n are301,406 and327,506. Within-layer sign/intersection proof is identical to the preceding audited affine-plane layer. Totals are therefore6,973,870,208 and12,397,783,154.

## Comparison

The independently read current repository README https://github.com/alexlegeartis/KissingNumbers on2026-11-17 displays87D6,226,194,610 and89D9,270,596,604. These are dated repository benchmarks, not an exhaustive worldwide literature-priority determination. The new89D count also supersedes the earlier independently audited ternary hyperplane count with its sparse layer11,622,942,176. The parent code and basic simplex embedding are classical; the finite head-assignment witnesses are the new work being checked.
