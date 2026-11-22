# Complementary-parity Steiner append in 91–93 dimensions

Add **40,960** rays in each of 91D, 92D and 93D. No earlier rays are removed.
The new bounds are **25,923,323,664**, **34,513,258,256**, and **47,398,160,144**.

The earlier triality and matching transfer uses Steiner quadruples of
22 rows. Each quadruple has a uniform axial mode; 348 also have both
uniform flat modes. On 20 of the remaining 37 supports, use a **mixed**
even/odd parity word and its complement. Each parity mode gets 64
zero-XOR axis patterns and all 16 sign patterns. Complementary modes
use different parity frames in every row. Distinct Steiner supports
overlap in at most two rows, so these 40 new modes do not interfere
with the existing Steiner stock. The mixed parity choices pass the
row-product bounds against all retained seven-row flat modes.
The appended points use no coordinates in the extra head.

`python3 verify.py` reads pinned earlier inventory and claims from
neighboring certificate folders, including the six-row Reed–Solomon
repack. The three columns of `witness.txt` are the Steiner support
index, the 22-bit support mask, and a four-bit odd-parity mask
indexed by increasing physical row.
