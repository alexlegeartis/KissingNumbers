"""Exact geometry and count transfer from the inherited checked sparse configurations."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys

if not __debug__:
    raise RuntimeError("Run without -O")
here = Path(__file__).resolve().parent
default = here.parent / "dim61-67-compatible-sign-families"
stock = Path(sys.argv[1]) if len(sys.argv) > 1 else default
for name, expected in {
    "verify.py": "5b7d63cb8038885889d825df67adb5af4596d81694c8876e057a7beff18a7286",
    "verification.json": "57e2ebbb6d626a10a9e9f876e19dce3d91ad83e6cb34801a70bbd194bd42b1ba",
}.items():
    assert sha256((stock / name).read_bytes()).hexdigest() == expected, name
old = json.loads((stock / "verification.json").read_text())
assert old["status"] == "PASS"
ancestor = {row["dimension"]: row["total"] for row in old["dimensions"]}

marks = (24, 42, 53, 55)
words = [sum(1 << j for j, x in enumerate(marks) if x & s == s)
         for s in range(64) if s.bit_count() <= 2]
span = {0}
for w in words:
    span |= {w ^ v for v in span}
assert span == set(range(16))

labels = (0, 1, 3, 7, 6, 4)
head = ((2, 0), (1, 1), (-1, 1), (-2, 0), (-1, -1), (1, -1))
for (i, (x, y)), (j, (u, v)) in product(zip(labels, head), repeat=2):
    assert 60 + x*x + 3*y*y == 64
    assert 60 - 2*(16-(i ^ j).bit_count()) + x*u + 3*y*v <= 32
    # Existing weight-4 vectors may use either head coordinate.
    for a, b in product(range(2), repeat=2):
        rem = 32 - 4*(4-a-b+a*abs(x))
        assert rem >= 0 and rem**2 >= 3*(4*b*abs(y))**2
    assert 8*abs(x) <= 32 and 3*(8*abs(y))**2 <= 32**2
assert 16*2 == 32

n = len((here / "witness.txt").read_text().split())
assert n == 80
dense = n << 22
assert ancestor[62] == (6 << 24) + 10610*2048 + 9145*16 + 124
assert ancestor[64] == (1 << 28) + 32768*2048 + 10416*16 + 128
assert ancestor[66] == (1 << 28) + (34028+1938)*2048 + 11088*16 + 132
assert ancestor[67] == 355090822 + 4624
modern = {d: ancestor[d] - (1 << 28) + dense for d in (64,66,67)}
modern[65] = modern[64] + 2
modern[62] = ancestor[62] - (6 << 24) + 6*(n << 18)
print(json.dumps({"status": "PASS", "marks": marks, "retained_labels": labels,
                  "rank": 4, "head_pair_types": 36, "dense64": dense,
                  "lower_bounds": dict(sorted(modern.items()))}, indent=2))
