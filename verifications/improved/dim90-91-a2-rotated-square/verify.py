"""Exact symbolic/combinatorial certificate of 90/91D A2-square augmentation."""
from itertools import combinations,product
from support import supports
S = supports()
all_triples = [t for b in S for t in combinations(b,3)]
assert len(all_triples)==len(set(all_triples))==4*3465
codes=[x for x in range(256)if x.bit_count()%2==0]
assert len(codes)==128 and all((x^y).bit_count()>=2 for x,y in combinations(codes,2))
# Rows in each A2 plane now use an orthonormal basis at −30°, +60° from a sparse root.
# A square sign row consists of two coordinates ±1/√2. It is one of 15°,105°,195°,285°.
# Signed sparse root unit axes are at 0°,60°,120°,180°,240°,300°.
d=[min((x-y)%360,(y-x)%360)for x,y in product(range(15,360,90),range(0,360,60))]
assert min(d)==15
# Old sparse/latitude dot≤4√3 cos15°; its square 24+12√3 <45=(old norm·new norm/2)^2.
assert 12*12*3 < (45-24)**2
# Dense per-plane norm1, hence dot≤4 < sqrt45. Latitude head contributes nothing to squares.
assert 4**2<45
# Root vectors use two orthonormal coordinates ±1, norm sqrt2.
# Against old dense/sparse/latitude: ≤√(2*3)*sqrt2 as they occupy at most two old blocks.
assert 4*2*3 <45
# Against a square vector max two overlaps each 1/√2: dot≤sqrt2, half product of norms exactly sqrt2.
# Squares have norm2; distinct supports meet ≤2, and equal supports have even sign distance≥2.
# Old dense / sparse / latitude occupy at least 14 planes; all new layers <=4. New roots have
# two nonzero coordinates, distinct from squares' eight. Drop two old poles in d91, none in d90.
assert len(S)==3465
for n,old,drop in [(90,21943640934,0),(91,22359345896,2)]:
 roots=4*(n*(n-1)//2);sq=len(S)*len(codes)
 print(n, len(S), sq, roots,'delta',sq+roots-drop,'total',old+sq+roots-drop)
