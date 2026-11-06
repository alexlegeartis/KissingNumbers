# Joint P48 coverage in dimensions 54–59

Six independent selections of existing isometric images improve the retained repository constructions by 1,066 points in total. These are explicit lower bounds, not optimality or worldwide-priority claims. The 60D result is unchanged.

| Dimension | Previous | New | Gain | Selected lines |
|---|---:|---:|---:|---:|
| 54 | 52,974,438 | 52,974,448 | 10 | 279,188 |
| 55 | 53,389,476 | 53,389,510 | 34 | 486,692 |
| 56 | 54,253,958 | 54,254,038 | 80 | 918,899 |
| 57 | 54,747,156 | 54,747,366 | 210 | 1,165,530 |
| 58 | 56,237,690 | 56,237,990 | 300 | 1,910,740 |
| 59 | 56,907,532 | 56,907,964 | 432 | 2,245,680 |

## Replay

Run `python verify.py` from any working directory using this script's absolute path. NumPy is required. The replay is read-only, takes about a minute on one CPU, and needs no network, GPU, search pool or outside coordinate file. `python final.py` derives the six counts from retained masks. `python test_certificate.py` checks rejection of altered certificates.

The verifier checks hashes, all 30,283,653 base pairs, every selected isometry, sign-canonical full-coordinate disjointness, and that ownership masks reproduce the full image union exactly. Signed E6/E7/E8 roots are generated explicitly. The 306-, 510- and 604-point tails are replayed in their exact integer/quadratic metrics, including every pair and antipodality. Each tail also supplies the poles.

## Construction and arithmetic

Take the 52,416,000-point unit P48 shell. Each selected mother line is owned by exactly one class, whose distinct signed heads have inner product at most 1/3. Assign one distinct antipodal tail pair to each class. Replace every selected signed head x by (sqrt(3)/2*x, +u/2) and (sqrt(3)/2*x, -u/2), and append the unit tail as poles. Same-class, different-class, same-head and pole products are bounded respectively by 1/4+1/4, 3/8+1/8, 3/4-1/4 and 1/2. Equator–cap products are at most sqrt(3)/4. Thus N = 52,416,000 + 2S + number_of_poles.

The only inherited geometric inputs are the published P48 minimum six and shell population 52,416,000, and the identification of the supplied Gram with that lattice. The finite tail and head claims are replayed here. No expanded all-pairs check of 50 million points is claimed.

## Provenance and reproducibility

Inputs are copied byte-for-byte from repository commit f0083954ce63e292b0d6024d73e35cb2486716e2: family_refilled.npz from dim60-antipodal840-lift/source, and capdirs_9/10/11.npz from dim49-63-p48-caps/data. manifest.json pins every input and certificate by SHA-256. The base contains 7,783 lines. Discovery selected among its existing 577 isometric images by positive single-image exchanges; results.json records those histories. Verification does not trust the discovery membership cache or its objective. It reconstructs the full images directly from the saved transforms, and independently derives first-owner masks and populations. Search-local optimality is not a prerequisite or a global maximality claim. No manuscript update is included.
