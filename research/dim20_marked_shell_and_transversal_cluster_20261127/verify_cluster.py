from collections import Counter
from itertools import combinations,product
# F4=F2[a]/(a²+a+1); 3=a+1.
def mul(a,b):
    p=0
    while b:
        if b&1:p^=a
        a<<=1
        if a&4:a^=7
        b>>=1
    return p
def add(a,b):return tuple(x^y for x,y in zip(a,b))
def times(t,a):return tuple(mul(t,x) for x in a)
def pv(a):return sum(x<<(2*i) for i,x in enumerate(a))
def dot(a,b):
    return mul(a[0],b[0])^mul(a[1],b[1])^mul(a[2],b[2])
X=list(product(range(4),repeat=3))
C=[(1,t,mul(t,t)) for t in range(4)]+[(0,0,1)]
U=(0,1,0)
def make_line(x,c):return tuple(sorted(add(x,times(t,c)) for t in range(4)))
def setmask(a):return sum(1<<pv(x) for x in a)
choices=[]
for i, c in enumerate(C):
    lines=sorted(set(make_line(x,c) for x in X));assert len(lines)==16
    tangent={add(times(t,c),times(s,U)) for t,s in product(range(4),repeat=2)}
    assert len(tangent)==16
    assert all(d not in tangent for j,d in enumerate(C) if i!=j)
    choices.append([(setmask(L)|setmask(M),(L,M)) for L,M in combinations(lines,2) if add(L[0],M[0]) in tangent])
    assert len(choices[-1])==24
def go(i,used):
    if i==5:return []
    for mask,row in choices[i]:
        if mask&used:continue
        rest=go(i+1,used|mask)
        if rest is not None:return [row]+rest
    return None
win=go(0,0); assert win
print('affine lines reps',[[L[0] for L in v] for v in win])
Bank=[]
for i in range(5):
    js=[k for k in range(5) if k!=i]
    needed=set(win[i][0]+win[i][1]); bank=[]
    for x in product(range(4),repeat=4):
        s=(0,0,0)
        for j,y in zip(js,x):s=add(s,times(y,C[j]))
        if s in needed:
            w=list(x);w.insert(i,-1)
            bank.append(tuple(w))
    assert len(bank)==32
    Bank.extend(bank)
for a,b in combinations(Bank,2):
    same=sum(x==y and x>=0 for x,y in zip(a,b))
    assert same<=2
print('PASS',len(Bank),'axis rows; all 16 signings each give',16*len(Bank),'same-profile transversal points')
