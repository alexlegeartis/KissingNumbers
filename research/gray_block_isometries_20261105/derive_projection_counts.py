"""Exact conditional section counts, not a proof of cross-layer compatibility."""
from fractions import Fraction
known={86:6225165512,87:6226194610,88:6921251584,89:9229354136,90:12306670328,91:16409599116,
       92:21880490288,93:31040958640,94:47697508076,95:77565848332,
       96:142162362560}
for d in range(86,97):
    if d<=92:
        m=96-d;n=2*Fraction(2**36*6**m,8**m)
    elif d<=94:
        m=97-d;n=2*Fraction(2**36*(6**m+2),8**m)+Fraction(2**37,2**m)
    elif d==95:
        n=2*Fraction(2**36*38,64)+2**35
    else:n=Fraction(3*2**36)
    assert n.denominator==1 and n>known[d]
    print(d,int(n),known[d])
