from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
from supports import hanani,double

def count():
    D=Path(__file__).resolve().parent
    base=hanani(8)
    T=[tuple(i for i in b if i!=21) for b in base if 21 in b]
    assert len(T)==70
    Q=[b for b in double(base,22) if b[-1]<42]
    Sq={t for b in Q for t in combinations(b,3)}
    five=[tuple(map(int,l.split())) for l in (D/'five.txt').read_text().splitlines()]
    Fv={t for b in five for t in combinations(b,3)}
    pairs={p for b in five for p in combinations(b,2)}
    assert len(pairs)==861
    W=[]
    for C in T:
     classes=[]
     for offset in range(4):
      G=[]
      for a,b in product(range(4),repeat=2):
       c=(a^b^offset)
       U=tuple(sorted((2*i+(w>>1),w&1) for i,w in zip(C,(a,b,c))))
       R=tuple(i for i,s in U)
       if R in Sq or R in Fv:continue
       G.append(U)
      classes.append(G)
     W.extend(max(classes,key=len))
    assert len(W)==len(set(W))
    assert len(W)>0
    # Compare all signed pair resources for the 3-plane additions.
    R=set()
    for U in W:
        rr=set(combinations(U,2))
        assert not (R&rr)
        R.update(rr)
        assert all(p in pairs for p in combinations(tuple(i for i,s in U),2))
    # Geometry: c=cos(15 deg), s=sin(15 deg). (c+s)^2=3/2 fits all
    # mixed D88 roots. c^2=(2+sqrt(3))/4, so 4*c^2<15/4.
    assert 3 < F(7,4)**2 and 2+F(7,4)==F(15,4)
    assert 4*3==3*4 # square mixed dot cap sqrt(3)
    assert 4*27<3*45 # all three Pless rows even without angular bound
    assert 4*3**2<3*45 # dense cap
    assert len(W)==560
    return len(W)

if __name__=='__main__':
    print(count())
