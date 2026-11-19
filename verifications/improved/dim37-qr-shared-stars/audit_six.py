"""Independent explicit verification using exact integer/surd products."""
from itertools import product,combinations
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'dim37-38-qr-merge'))
from supports37 import supports

def le_surd(a,b,c):
    """Check a+b sqrt(2)<=c in the nonnegative b cases used here."""
    assert b>=0
    return a<=c and b*b*2<=(c-a)**2
R=json.loads((P/'six_recipe.json').read_text())['ordinary_supports']
d=json.loads((P/'star_recipe.json').read_text())
Z=[(36,tuple(x)) for x in d['partition']]+[(int(i),tuple(x)+(36,)) for i,x in d['ordinary'].items()]
W=[set(i for i in range(36) if w>>i&1) for w in R]
A=supports()
assert all((x&a).bit_count()<=4 for x in R for a in A)
for h in W:
    assert len(h)==6
    for c,L in Z:
        a=8*(c in h)+4*len(h&set(L));b=8 if c==36 else 4
        assert le_surd(a,b,16)
assert {i for c,L in Z if c==36 for i in L}==set(range(36))
# Dense cross squared: (24+2m sqrt(2))^2 = 576+8m^2+96m sqrt(2).
for m in (1,4):
    assert le_surd(576+8*m*m,96*m,8*(144+m*m))
# Explicitly view each row as (positive-bit mask, support mask, signed radical).
pts=[]
for w,h in zip(R,W):
    h=sorted(h)
    for sign in product((-1,1),repeat=6):
        if sign.count(-1)&1:continue
        m=sum(1<<i for i,s in zip(h,sign) if s<0)
        pts.extend([(w,m,1),(w,m,-1)])
assert len(pts)==64*len(set(R))==len(set(pts))
def product0(x,y):
    a,b,u=x;c,d,v=y;k=a&c
    return 4*(k.bit_count()-2*((b^d)&k).bit_count())+8*u*v
assert max(product0(x,y) for x,y in combinations(pts,2))<=16
report=dict(status='PASS',dimension=37,total=534008+len(pts),added_seven_site=len(pts),
            seven_seven_pairs=len(pts)*(len(pts)-1)//2)
print(json.dumps(report,sort_keys=True))
