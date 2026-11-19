"""Independent sparse sign and signed-root accounting in 37D."""
from collections import Counter
from itertools import combinations,product
from pathlib import Path
import sys,json
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'dim37-38-qr-merge'))
from supports37 import supports
d=json.loads((P/'star_recipe.json').read_text())
sets=[(36,tuple(b)) for b in d['partition']]+[(int(a),tuple(b)+(36,)) for a,b in d['ordinary'].items()]
assert len(sets)==len(set(sets))==45
# Each nonzero leaf magnitude is 2. Each centre magnitude is 4.
R=supports()
block=[]
for c,l in sets:
    assert len(l)==len(set(l))==4 and c not in l
    assert all(0<=i<37 for i in (c,)+l)
    r=sum(1<<i for i in l)
    assert all(8*(s>>c&1)+4*(r&s).bit_count()<=16 for s in R)
    for a in (1,4):
        w4=a if c==36 else 2
        q=4*w4+sum(2*(a if i==36 else 2) for i in l)
        assert 4*q*q<=32*(144+a*a)
    for sign in product((-1,1),repeat=4):
        if sign.count(-1)&1:continue
        for u in (-1,1):
            block.append(dict([(c,4*u)]+[(i,2*t) for i,t in zip(l,sign)]))
assert len(block)==len({tuple(sorted(z.items())) for z in block})==720
m=max(sum(v*y.get(i,0) for i,v in x.items()) for x,y in combinations(block,2))
assert m<=16
def edge(c,l):
    return (min(c,l),max(c,l))
def incidences(S):
    return {edge(c,j) for c,v in S for j in v}
# If the root is not incident to a centre and one of its leaves, its
# two magnitudes on any star sum to at most four. Each possible sign
# pair is blocked for every incident unsigned edge because all centre
# signs and either leaf sign occur over the even sign words.
I=incidences(sets)
assert len(I)==144
P0=incidences([(c,l) for c,l in sets if c==36])
assert len(P0)==36 and len(I-P0)==108
n=180224+128*len(R)+2*37*36+len(block)-4*len(I)
assert n==534008
print(json.dumps({'status':'PASS','dimension':37,'total':n,'signed_roots_removed':4*len(I),
 'star_pairs':len(block)*(len(block)-1)//2,'max_star_star':m},sort_keys=True))
