"""Direct exact generator and distance of binary DG(6,1) length64."""
from itertools import combinations,product
from functools import reduce
from operator import xor
def code():
    def mul(a,b):
        t=0
        for _ in range(5):
            if b&1:t^=a
            b>>=1;a<<=1
            if a&32:a^=37
        return t
    assert all({mul(a,b) for b in range(32)}==set(range(32)) for a in range(1,32))
    def tr(a):
        t=0
        for _ in range(5):t^=a;a=mul(a,a)
        assert t<=1
        return t
    B=[1<<i for i in range(5)]
    F=[]
    pairs=list(combinations(range(6),2))
    for a,b in product(range(32),repeat=2):
        def Q(x,y):return tr(mul(a,mul(x,y))^mul(b,mul(mul(x,x),y)^mul(x,mul(y,y))))
        d=[Q(x,x) for x in B]
        h=sum((Q(B[i],B[j])^d[i]*d[j] if j<5 else d[i])<<r for r,(i,j) in enumerate(pairs))
        F.append(h)
    assert len(set(F))==len(F)==1024
    D={a^b for a,b in combinations(F,2)}
    assert len(D)==15903
    row=lambda t:sum(1<<b for b in range(64) if (b&t)==t)
    E=[row((1<<i)|(1<<j)) for i,j in pairs]
    R=[row(1<<i)for i in range(6)]+[2**64-1]
    val=lambda d:reduce(xor,(x for i,x in enumerate(E) if d>>i&1),0)
    A=[0]
    for v in R:A += [x^v for x in A]
    assert len(set(A))==128 and min(x.bit_count() for x in A if x)==32
    min_d=64
    for u in D:
        z=val(u)
        for y in A:min_d=min(min_d,(y^z).bit_count())
    assert min_d==24
    return [val(x) for x in F],A
if __name__=='__main__':
    G,A=code()
    print(len(G),len(A),len(set(G))*len(set(A)))
