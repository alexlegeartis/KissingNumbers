from math import comb
from fractions import Fraction as F
M={1:F(5517,340),2:124,3:F(83669,170),4:899,9:F(35061,68)}
def C(a,b):return comb(a,b) if 0<=b<=a else 0
def K(j,d):
    return sum((-1)**a*C(d,a)*C(32-d,j-a) for a in range(j+1))
def s(d):
    return 1+sum(F(w)*F(K(j,d),comb(32,j))for j,w in M.items())
assert s(0)==2048 and min(M.values())>=0
assert all(s(d)<=0 for d in [12,16,20,24,32])
print("Exact Krawtchouk bound",s(0))
