"""Independent Cartesian enumeration of every relevant D4-block tail vector."""
import verify as p
from collections import Counter
from itertools import product
import numpy as np
from sympy import Matrix

B=Matrix([[1,-1,0,0],[0,1,-1,0],[0,0,1,-1],[0,0,1,1]])
assert B*B.T==Matrix(p.d4.tolist()) and abs(B.det())==2
J=np.array((2*B.inv()).tolist(),dtype='i8')
def cartesian(r):
    out=[]
    for row in product(range(-2,3),repeat=4):
        q=sum(t*t for t in row)
        if q>8 or sum(row)%2:continue
        z=np.array(row,dtype='i8')@J
        assert (z%2==0).all()
        z=z//2
        if z[r:].any():continue
        out.append((q,z[:r]))
    return out

for r in range(1,11):
    blocks,raw=p.candidate(r)
    polys=[Counter({(0,(0,)*72):1})]
    offset=0
    for v in blocks:
        s=len(v);X=raw[offset:offset+s,:];offset+=s
        D=Counter()
        for q,z in cartesian(s):
            D[q,tuple(z@X%2)]+=1
        new=Counter()
        for (a,d),k in polys[-1].items():
            for (b,e),l in D.items():
                if a+b<=8:new[a+b,tuple(x^y for x,y in zip(d,e))]+=k*l
        polys.append(new)
    C=Counter()
    for (q,code),c in polys[-1].items():
        if q:
            C[q]+=c*p.contact(code,8-q)
    own=p.report[r-1]
    assert all(C[q]==own['contacts_by_vertical_norm'][q]for q in (2,4,6,8))
    assert 6218175600+sum(C.values())==own['count']
    print('CARTESIAN PASS',72+r,dict(C))
