#!/usr/bin/env python3
"""Replay the QR48 fibres and exact surd gates. Run: python verify.py."""
import json
from collections import Counter
p=47
R=sum(1<<i for i in {j*j%p for j in range(1,p)})
mask=(1<<p)-1
basis={}
for j in range(p):
    w=((R<<j) | (R>>(p-j))) & mask
    w |= (w.bit_count() & 1)<<47
    while w:
        k=w.bit_length()-1
        if k not in basis:
            basis[k]=w
            break
        w ^= basis[k]
assert len(basis)==24
G=list(basis.values())
x=0
counts=Counter({0:1})
prefix=Counter({0:1})
dcentral=40
minimum_dense_difference=40
for i in range(1,1<<24):
    x ^= G[(i & -i).bit_length()-1]
    w=x.bit_count()
    counts[w]+=1
    h=x & 255
    prefix[h]+=1
    tail=w-h.bit_count()
    if h.bit_count()<=1: dcentral=min(dcentral,tail)
    if h.bit_count()<=2: minimum_dense_difference=min(minimum_dense_difference,tail)
assert min(k for k in counts if k)==12
assert len(prefix)==256 and set(prefix.values())=={65536}
assert dcentral>=11 and minimum_dense_difference>=10
assert sum(prefix[0 if j==8 else 1<<j] for j in range(9))==589824
# t=a+1=sqrt(41)-1, a=sqrt(41)-2; D=2(a+1)^2.
# Surd A+B*sqrt(41) -> exact positivity.
def ge(A,B):
    if A>=0 and B>=0:return True
    if A<0 and B<=0:return False
    return A*A>=41*B*B if A>=0 else 41*B*B>=A*A
# d >= 11 gives donor-dense dot <= a+19.
# 10D-(a+19)^2 = 362-74a = 510-74*sqrt(41)>0.
assert ge(510,-74)
# a>=3 and hence max(8,a+5)<=2(a+1)
assert ge(-5,1)
new=589824-1+40+4267*128+2*40*39
assert new==1139159
print(json.dumps(dict(qr48_weight_distribution=dict(sorted(counts.items())),
    first8_prefix_rank=8,central_to_retained_dense_Hamming=dcentral,
    projected_minimum_for_all_9_prefix_fibers=minimum_dense_difference,
    dense_original=589824,parents_removed=1,daughters=40,
    signed_octads=4267*128,roots=3120,total=new),indent=2))
