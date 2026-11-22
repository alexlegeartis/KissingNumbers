"""Exact two-site hyperplane section of the 80-cubic-coset 64D construction."""
import hashlib
from pathlib import Path
import json
import sys
from itertools import combinations
import numpy as np

if not __debug__:
    raise RuntimeError("Use Python without -O")
here = Path(__file__).resolve().parent
parents = Path(sys.argv[1]) if len(sys.argv) > 1 else here.parent
dense = Path(sys.argv[2]) if len(sys.argv) > 2 else here
src = parents / "dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz"
cert = parents / "dim61-67-compatible-sign-families/data/global64.json"
pins = {src: "cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6",
        cert: "8c0c35f8127c578118f1f09e5950d712a272db32f2d50f382e6fe32f46c6dbf5",
        dense / "witness.txt": "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f"}
for p, expected in pins.items():
    assert hashlib.sha256(p.read_bytes()).hexdigest() == expected, p
with np.load(src,allow_pickle=False) as z:
    old = list(map(int,z['words']))
payload=json.loads(cert.read_text())
new=[u["support"] for u in payload["new"]]
rows=old+new
assert (len(old),len(new),len(set(rows))) == (30828,1940,32768)
assert all(0 < s < 1 << 64 and s.bit_count() == 16 for s in rows)
i,j=2,41
m=(1 << i)|(1 << j)
none=sum(s&m==0 for s in rows)
both=sum(s&m==m for s in rows)
assert (none,both)==(18352,1960)
# All [16,11,4] classes, modified or new, have a dual of min weight 8.
# Hence selecting the two equal bits keeps half of a class with both sites occupied.
classes = {u["support"]: u["checks"] for u in payload["old"] + payload["new"]}
assert len(classes) == len(payload["old"]) + len(payload["new"])
stock_checks = [65535] + [sum(1 << z for z in range(16) if z >> b & 1) for b in range(4)]
def chart(checks):
    span = {0}
    for h in checks:
        span |= {h ^ x for x in span}
    assert len(span) == 32 and {x.bit_count() for x in span} == {0,8,16}
chart(stock_checks)
for s in rows:
    if s&m == m and s in classes: chart(classes[s])
sparse=none*2048+both*1024

# Direct affine stereotyping for the four-support words (all 16 signs).
blocks={sum(1<<x for x in (a,b,c,a^b^c)) for a,b,c in combinations(range(64),3)}
assert len(blocks)==10416
quad=sum(16 if b&m==0 else 8 if b&m==m else 0 for b in blocks)
poles=2*(64-2)
q3=list(map(int,(dense/"witness.txt").read_text().split()))
flags={0}
for s in range(64):
    if s.bit_count() > 2: continue
    x=int(i&s==s) + 2*int(j&s==s)
    flags|={x^u for u in flags}
assert flags==set(range(4)) and len(q3)==80
dense63=len(q3)*(1<<21)
assert (sparse,quad,poles,dense63)==(39591936,146568,124,167772160)
print(json.dumps(dict(status="PASS",dimension=63,coordinate_pair=[i,j],dense=dense63,
    sparse=sparse,quadruple=quad,poles=poles,total=dense63+sparse+quad+poles),indent=2))
