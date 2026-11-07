# Individual base-triangle switches in 36D and 37D

This frozen exact construction gives **501,813 points in 36D** and **505,629 in 37D**, an increase of 126 in each dimension over the retained completion. It switches 90 three-point base triangles, adds 1,511 equatorial phase orbits, removes 1,469, and retains all 636 latitude body orbits at heights ±2. The base cardinality remains 403,704.

`proposal.json` is the immutable discovered delta. `switches.json` independently records each switched large vector and its source/target coordinate blocks. `certificate.json` pins hashes, counts, and proof scope. `proof.md` describes the exact construction.

The two verifier routes are:

- `verify.py`: fresh enumeration using the retained standard-library rational shell enumerator, complete base reconstruction, exact compressed block-product tables, every affected pair category, distinctness, and exact ranks 36/37.
- `verify_coordinates.py`: an independent rational LDL enumeration with no producer/search/source-enumerator imports, complete base reconstruction, bounded-int64 direct products for every affected pair category, and full 36D/37D norms and uniqueness.

Each standalone verifier freshly replays the unchanged completion proof in isolation. `negative_controls.py` restores deleted source points and confirms products 8 and 10, above the allowed 6. The shipped receipts are reproducible checks, not assumptions. The root package wrapper runs every script in temporary copies so reports do not invalidate the manifest.

Run `python ../verify.py` for the complete immutable package replay. For optional explicit coordinates, run `python verify.py --write-coordinates --prefix exported` in a writable copy. The 36D coefficient rows use the supplied positive-definite Gram matrix diag(G,G,G), squared norm 12 and maximum distinct product 6. The 37D file provides 1,908 norm-eight body vectors and both heights ±2; normalize all full points by sqrt(12).

Source `../completion/certificate.json` remains byte-for-byte unchanged. Its count is 501,687/505,503; older `../source/vectors.npz` contains an earlier equator and is not used as the current phase code. All source short-shell populations are enumerated exactly. No large-shell completeness assumption or numerical angle tolerance is added. The paper is unchanged, and worldwide priority is not asserted beyond the dated source comparison.
