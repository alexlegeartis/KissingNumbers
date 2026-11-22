# Triality and matching-row transfer in dimensions 91 to 93

| Dimension | Predecessor | New lower bound | Gain |
|---:|---:|---:|---:|
| 91 | 25,917,421,568 | 25,922,303,504 | 4,881,936 |
| 92 | 34,507,356,160 | 34,512,238,096 | 4,881,936 |
| 93 | 47,392,258,048 | 47,397,139,984 | 4,881,936 |

This modifies only the flat augmentation in the published [six- and seven-row configurations](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim91-93-flat-parity-extensions). All 24/32/44 individually normalized Z8 fibers and all 145,004,544 Golay sparse points stay unchanged. Among the predecessor's 2,613,248 six- and seven-row flats, retain exactly 131,072 points, represented by 16 seven-row parity modes. Add 1,107,472 triality-frame and ordinary D4 points and 6,256,640 matching-frame points. These make a 7,495,184-point replacement for the former flat stock.

The ordinary body is identical to the 89D Z8/Golay body: 22 four-coordinate rows, each dense row in one of the even-flat, odd-flat or axial 24-cell frames, and body squared radius 88. Append every transferred point at zero in the new head, so each dense squared radius is at least 88 even after adding its own nonzero head. Different ordinary occupancy or row shapes separate the layers. All finite data are decoded using the exact masks and frame conventions from the [published 89D matching witness](https://github.com/alexlegeartis/KissingNumbers/tree/main/verifications/improved/dim89-mixed-matching-rows). The supplementary 352-word length-eight inner axis code is due to [Laaksonen and Östergård](https://www.cs.helsinki.fi/u/ahslaaks/codes/), and its minimum distance is replayed exactly.

Use the 385 quadruple supports from the Steiner quadruple system of order 22. Each support carries the uniform axial mode, and selection.json gives 696 additional frame-parity modes with a uniform even or odd flat frame. Each fixed mode uses four orthogonal projective axes per row, a 64-word zero-XOR axis code and all 16 signs, for 1,024 rays at squared radius 16. Add the 24 D4 roots in each ordinary row. The transferred matching patterns have 2/850/3/61 supports of weights 5/6/7/8, respectively, with the same inner axes and signs as in 89D.

For the old flat-mixed check, let q be their overlapping rows. Between a new matching pattern and any old flat pattern, the raw dot is at most 2√2 q. A matching pattern of weight w and an old flat of weight r have squared radii 4w and 4r, so the checker verifies 2q²≤wr. For a new triality frame and an old flat parity mode, the raw cap is 4q−2d, where d is the number of distinct parity frames on common rows. This gives the exact integer test (4q−2d)²≤16r. Any uniform axial mode has mixed raw dot at most 2q≤8, and each D4 root has mixed raw dot at most 2. The other pair checks, including new matching–matching overlap plus equal tags, Golay comparisons, and disjointness from dense fibers, are repeated by the verifier or inherited unchanged from the cited parents.

Run python3 verify.py with assertions enabled. The original flat support ledger, predecessor claims and Golay bootstrap files are hash-bound to the inherited witnesses. Reproduction does not need the optimizer trade.py; it only uses the explicitly stored selection.
