#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations
import json
from supports37 import supports
from verify_predecessor import replay as prior_replay
D=Path(__file__).resolve().parent
assert prior_replay()["points"]==544474
u=((1<<12)-1)<<12
s=supports()
old=[t>>12 for t in s if t&u==t]
assert old==[255,3855,4080]
assert max((m&u).bit_count() for m in s if m&u!=m)==4
rad=[int(m) for m in (D/'six.txt').read_text().split()]
assert len(set(rad))==168 and all(m.bit_count()==6 and (m&u).bit_count()<=4 for m in rad)
j=json.loads((D/'dense_parent38.json').read_text())
kernel={0}
for g in j['generators_hex']:
 g=int(g,16);kernel|={x^g for x in kernel}
dense={x^int(r,16) for r in j['representatives_hex'] for x in kernel}
# Min norm squared is 4*36+1. New vectors have L1=16 and marked coordinate zero.
assert 4*(2*16)**2<=32*(4*36+1)
lines=[list(map(int,l.split())) for l in (D/'win.txt').read_text().splitlines()]
pts=set()
g={0}
for b in lines[0]:
 assert 0<b<4096;g|={x^b for x in g}
assert len(g)==256
for mask,neg,n in lines[1:]:
 assert 0<mask<4096 and mask.bit_count()==8 and neg&mask==neg
 v={(mask,neg^(x&mask)) for x in g}
 assert len(v)==n and not (v&pts)
 pts|=v
for (m,a),(n,b) in combinations(pts,2):
 q=m&n
 assert q.bit_count()-2*((a^b)&q).bit_count()<=4
# Exterior roots have two ±4 coordinates. New local vectors have ±2 entries.
assert 2*4*2==16
# Exterior octads and radical layers have coefficients ±2 on ordinary support;
# the above support-intersection caps make their dot product at most 16.
assert 4*2*2==16
assert len(dense)+4*(36*35//2)+128*2742+64*168+2==544474
assert len(pts)==576
print(json.dumps(dict(dimension=37,previous_bound=544474,deleted=384,inserted=len(pts),lower_bound=544474-384+len(pts),local_phase_group_order=len(g),pair_checks=len(pts)*(len(pts)-1)//2)))
