# Nonlinear sign-code refinement of the audited RS9 sparse layer

Retain exactly the support-ten extended Reed–Solomon construction and D_n roots proved and independently audited in research/rs9_sparse_20261117. Replace the64 binary Hamming sign words on every support by the72 explicit length10 binary words in certificate.json. These are the first code from the primary authors' electronic supplement https://users.aalto.fi/~pat/opt/72.txt to P. R. J. Östergård, T. Baicheva, and E. Kolev, “Optimal binary one-error-correcting codes of length10 have72 codewords,” IEEE Transactions on Information Theory45(1999),1229–1231; explanatory source https://users.aalto.fi/~pat/72.html .

The verifier checks every one of the2556 sign-word pairs and finds minimum Hamming distance3. Thus same-support inner products remain≤18. Distinct-support, root, and dense cross bounds are unchanged because they used only absolute coordinate values and support intersections, not sign symmetry.

IMPORTANT: the72-word sign code is not antipodal, and NO extra negatives are appended. The complete kissing configuration need not be antipodal. The dense code remains antipodal, and its interaction with these sparse points is covered by the existing absolute cross bound. The sparse count is exactly72*k*9^5+2n(n−1), not twice this number.

The new sign replacement is pending independent audit. All underlying RS9 geometry was already independently audited. This remains an additive improvement, not a large-factor breakthrough.
