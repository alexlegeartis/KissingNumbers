"""Exact regression: a partially overlapping translate union can remain valid."""
from pathlib import Path
import json,time
OUT=Path(__file__).resolve().parent
start=time.monotonic()
def gray(a):
 return sum((((a>>3)&1)^((a&1)*(x&1))^(((a>>1)&1)*((x>>1)&1))^(((a>>2)&1)*((x>>2)&1)))<<x for x in range(8))
T=[gray(a)for a in range(16)]
C={T[a]|T[b]<<8|T[c]<<16|T[(a+b+c)%16]<<24 for a in range(16)for b in range(16)for c in range(16)}
translated={v^0x0f0f0f0f for v in C};points=sorted(C|translated)
if (len(C),len(C&translated),len(points))!=(4096,2560,5632):raise ValueError('counterexample cardinalities changed')
minimum=32;checked=0
for i,x in enumerate(points):
 for y in points[i+1:]:
  d=(x^y).bit_count();checked+=1
  if d<8:raise ValueError('invalid distinct pair')
  minimum=min(minimum,d)
if minimum!=8 or checked!=15856896:raise ValueError('counterexample distance or coverage changed')
receipt=dict(status='PASS',scalar=1,repeated_local_mask=15,repeated_binary_mask_hex='0f0f0f0f',parent_words=len(C),intersection_words=len(C&translated),union_words=len(points),minimum_distinct_distance=minimum,distinct_pairs_checked=checked,elapsed_seconds=time.monotonic()-start,implication='minimum cross distance0 proves overlap, not an invalid union; unrestricted262144 bound retracted')
(OUT/'overlap_counterexample.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
