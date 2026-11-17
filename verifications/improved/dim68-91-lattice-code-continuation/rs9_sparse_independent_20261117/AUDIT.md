# Independent extended-F9 Reed–Solomon sparse layer audit

PASS, conditional only on the already-audited dense QR48 configurations:
87D6,996,258,582;88D9,324,561,000;89D12,427,704,400;90D20,954,734,650;91D20,954,735,010.

This sparse layer REPLACES the preceding affine-plane/affine-flat layers. The counts must not be added to those previous sparse totals. Dense counts are2*M*3^18 with M=9,12,16,27,27, and norm squared45. The published QR48 minimum premise remains as documented in the preceding dense-head audits; the sparse construction itself has a complete independent finite verification.

## Support construction and independent exhaustive check

Use nine disjoint groups of nine axes, indexed by t in F9. A tenth group has k=min(9,n-81) axes. For every degree<=5 polynomial overF9 with leading coefficient among the selected k symbols, use its evaluation label in each finite group and its leading coefficient in the tenth group. Each support has ten coordinates. There are k*9^5 supports. Different polynomials have at most five agreeing positions: if their leading coefficients agree the difference has degree<=4 and adds one infinity agreement; otherwise degree<=5 and infinity disagrees. Nine evaluations distinguish degree<=5 polynomials.

The auditor additionally enumerated ALL9^6=531,441 extended evaluation codewords over F9=F3[u]/(u^2+1), using independently implemented field tables and Horner evaluation. Associativity, distributivity and inverses were checked on all field elements. Every nonzero evaluation word has Hamming weight at least5, with minimum exactly5. Linearity then independently establishes support intersection<=5. This is much smaller than the QR48 enumeration and was fully executed, not inferred from a search timeout. The exact weight histogram is in report.json.

## Signs and internal geometry

The ten parity-check columns are the nonzero four-bit columns except{1,2,4,8,15}. Independent enumeration of all1024 binary vectors reproduces64 zero-syndrome words, and every pair has Hamming distance at least3. The sign code is closed under complement, so antipodes occur exactly once within its64 labels. Each support coordinate has magnitude sqrt(45/10). Same-support distinct signs have inner product<=18; different supports have inner product<=5*(45/10)=22.5. Distinct supports or sign vectors cannot duplicate.

Add the2n(n-1) scaled D_n roots, each with two entries+/-sqrt(45/2). Their mutual inner products are<=22.5. Root/support-ten cross inner product is at most45/sqrt5<22.5, independently checked by rational squaring. Support-ten and root vectors have different support sizes and are disjoint.

## Dense cross bounds

For compressed87–89D, the84 body coordinates precede all head coordinates. Thus all heads lie in the final group, and a ten-support contains at most one head coordinate. Body entries have absolute value<=1; head entries<=sqrt(3/2). For90D and91D the bound remains valid, with no nonzero head and an optional unused zero axis. A dense/support-ten scalar product is at most(9+sqrt(3/2))*sqrt(45/10)<22.5. Squaring and using the exact value( sqrt(3/2) )^2=3/2 reduces this to sqrt(3/2)<5/3, verified by3/2<25/9. Dense/root cross products are at mostsqrt135<22.5. Bounds are absolute, so they also cover antipodes.

Every dense vector has at least42 nonzero coordinates, while sparse support sizes are10 or2. Hence there are no cross-layer duplicates. The final91D support family uses only90 axes; D91 roots may use the extra axis and remain valid.

## Scope and artifacts

verify.py is an independently written exhaustive F9 evaluation/sign verifier. input.json freezes the producer certificate, report.json records all counts and its hash, and replay.log records execution. There is no claim of worldwide novelty. The field-code and sparse-support ingredients are classical; this audit validates their use with the newly audited dense configurations.
