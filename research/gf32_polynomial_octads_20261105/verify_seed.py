"""Verify the structural 32D seed used by the polynomial-support search."""
import json,itertools
from pathlib import Path
c=json.load(open(Path(__file__).resolve().parent/'seed32.json'));F=c['full_even_supports'];E=c['exceptional_antipodal_lines'];G=c['dense_generator'];d=32;mask=(1<<d)-1
assert len(F)==1704 and len(set(F))==len(F) and all(type(s) is int and s.bit_count()==8 and not(s>>d) for s in F)
assert all((s&t).bit_count()<=4 for s,t in itertools.combinations(F,2))
assert len(E)==816 and len(set(map(tuple,E)))==len(E)
for s,w in E:
 assert s.bit_count()==8 and not(s>>d) and not(w&~s) and not(w&(s&-s))
 assert all((s&t).bit_count()<=4 for t in F)
for (s,v),(t,w) in itertools.combinations(E,2):
 ov=s&t;assert abs(ov.bit_count()-2*(ov&(v^w)).bit_count())<=4
assert len(G)==17 and all(type(x) is int and 0<x<=mask for x in G)
span=[0]
for g in G:span += [x^g for x in span]
assert len(set(span))==1<<17 and min(x.bit_count() for x in span if x)>=(d//4)
# 128 even signings per full octad; 2 per exceptional line;
# 2^17 dense binary signs; 2*d*(d-1) two-coordinate roots.
N=128*len(F)+2*len(E)+len(span)+2*d*(d-1)
assert N==c['points']==352800
print('PASS',N,'points: source octad constraints, binary-code minimum distance, and disjoint layers')
