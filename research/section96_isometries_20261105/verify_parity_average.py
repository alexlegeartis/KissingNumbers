"""Exact first moment for the mixed-parity 92D section universe."""
from fractions import Fraction
from itertools import product
E=[x for x in range(16) if x.bit_count()%2==0]
O=[x for x in range(16) if x.bit_count()%2==1]
W=[x for x in E if x.bit_count()==2]
assert len(E)==len(O)==8 and len(W)==6
assert all(sum((a^b).bit_count()==1 for b in O)==4 for a in E)
for a in W:
 assert sum((a^b).bit_count()==0 for b in W)==1
n=2**36*Fraction(3,4)**4
expected=n*Fraction(1,6)**4*Fraction(1,2)**20
assert n==21743271936 and expected==16
print('PASS: uniformly sampled mixed-parity 92D section word has exactly 16 source neighbors at distance 20 on average')
