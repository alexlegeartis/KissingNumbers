# Matching-row append in dimensions 94 and 95

| Dimension | Prior lower bound | New lower bound | Gain |
|---:|---:|---:|---:|
| 94 | 74,138,716,876 | **74,139,180,748** | 463,872 |
| 95 | 117,088,390,212 | **117,088,854,084** | 463,872 |

The same append works on the 22 ordinary four-coordinate rows in
both dimensions, with zero head. It retains every point of the
previous [hexary parity configuration](../dim94-96-hexary-parity-packing/).
The numerical witness gives one four-row, 55 six-row, and 29
seven-row tagged patterns. Each tag fixes a pairing of the row
coordinates and its four orthogonal root axes, with two coordinates
of magnitude √2. The inner axis code and sign code give 1,024,
4,096, or 8,192 points per pattern respectively.

For shared rows, the dot bound depends on the number of overlaps
plus the number with matching tags. The checker compares all
new pairs and all 1,064 shifted Golay patterns. It replays the
BCH-flat and six-row parity support sets, and checks the ordinary
intersections with the flat and amplitude-four Reed–Solomon layers
and the full root layer. The 94D and 95D amplitude-four RS layers
induce the same coordinate-group partition on the ordinary body.
Each ordinary dense row gives dot at most 2√2 and both old dense
radii squared equal 96. Distinct row shapes, support sizes and
head coordinates separate the append from the preceding layers.

Run `python verify.py` with assertions enabled. This finite
certificate takes no optimizer or downloaded inputs.
