Dimension 28: 205,232 points (24 October 2026)

This explicit kissing configuration improves Kevin Li's 205,172-point
construction by 60. It is a modest new lower bound, not a large breakthrough
across dimensions 25--31. No improved bound in the other six dimensions is
claimed by this artifact.

Credit and provenance
--------------------
The Leech lifting method and the 512 odd heads in each of eight blocks are
from Kevin Li's public construction:
https://github.com/Felpix-Studios/kissing-number/tree/main/preprints/leech-lift-d28-d31
The full source file SHA-256 is
7e81d26b83b912549e4ba078691cd74f67350d51794edd6f934d210eb97a2734.
The downloaded repository archive SHA-256 is
eda4bab96d7d68be1629e7e18f15b5fd19d4ebf2d7986ff6d8fa8bf6fc54de56.
Source data and the unmodified all-pairs checker are provided under the
Apache 2.0 license; see THIRD_PARTY_LICENSE.txt. The new search jointly
chooses the even heads and the auxiliary layer. This discovery was made
with Codex in the KissingNumber27 research workspace.

What changed
------------
Four head sets have 524 vectors and four retain 528, for 4,208 owners.
The auxiliary layer grows from 164 to 256 points. Returning 16 owners to
the equator and removing their three caps costs 32 net points. The 92
additional auxiliary points therefore give a net gain of 60.

             old             new
equator      192,336         192,352
caps          12,672          12,624
auxiliary        164             256
total        205,172         205,232

All coordinates are exact integer combinations of 1, sqrt(2), sqrt(3),
and sqrt(6). Each row represents its normalization to the unit sphere.
The provenance JSON records the selected head vectors and tail sums.

Reproduce the proof
-------------------
Python 3.9+ and NumPy suffice. No download or original workspace is needed:

  python verify.py D28-205232.txt

This independently written checker rebuilds the Golay Leech shell, checks
the complete partition and all head-set and tail conditions, and checks
every interaction involving the auxiliary layer with unbounded-integer
arithmetic. It uses the standard kissing property of the Leech minimal
shell for equator/equator pairs. No tolerance is used in an acceptance
decision. The expected final status is EXACT VERIFICATION PASSED.

For an independent, exhaustive all-pairs check with Kevin Li's checker:

  g++ -O2 -std=c++17 -pthread check_kissing.cpp -o check_kissing
  ./check_kissing D28-205232.txt --threads 6

The C++ checker checks every pair, including all equator/equator pairs.
Its floating-point maximum cosine is informational only; acceptance uses
exact arithmetic. The result, when present, is recorded in
all_pairs28_205232.log. Compiler/platform timing varies substantially.

Scope of the search
-------------------
The search fixes the odd heads and varies compatible even Leech owners
and a finite pool of norm-eight tetrad/coordinate-axis heads over dual-D4
tail directions. Its best feasible objective is 480; its latest solver
upper bound is 530 for this restricted finite model. Neither is an upper
bound for the unrestricted kissing number. We do not claim optimality of
the new configuration or of this construction family.
