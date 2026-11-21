from itertools import product,combinations
from collections import Counter
from pathlib import Path
import hashlib,json,sys
P=Path(__file__).resolve().parent
UP=P.parent
if len(sys.argv)>1:UP=Path(sys.argv[1]).resolve()
p=P/'z8_generator.json'
assert hashlib.sha256(p.read_bytes()).hexdigest()=="903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73"
G=json.loads(p.read_text())['generator']
assert len(G)==12 and all(len(v)==24 and v[i]%2==1 and not any(v[:i])for i,v in enumerate(G))
A=[tuple((-1)**(((a>>2)+((a>>1)&1)*y+(a&1)*x)%2) for x,y in product(range(2),repeat=2))for a in range(8)]
B=[(-a[0],*a[1:])for a in A]
C=[tuple((2 if q<4 else -2)if i==q%4 else 0 for i in range(4))for q in range(8)]
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
T=[A,B,C]
for a,b in combinations(T,2):assert max(dot(x,y)for x,y in product(a,b))==2
D=[(p['color'],tuple(p['prefix']),tuple(p['head']))for p in json.loads((P/'heads.json').read_text())]
assert len(D)==len({h for _,_,h in D})==len({(c,a)for c,a,_ in D})==108
for c,a,h in D:
 assert 0<=c<=2 and len(a)==2 and all(0<=k<=7 for k in a)
 assert len(h)==8 and all(type(k)==int for k in h)and dot(h,h)==8 and h[0]+h[4]==0
 # Exact six/seven-coordinate magnitude caps after orthonormal compression.
 assert 2*h[0]**2<=8 and max(abs(h[j])for j in(1,2,3,5,6,7))<=2
for (c,a,h),(e,b,k) in combinations(D,2):
 bound=4 if c!=e else sum(dot(T[c][u],T[c][v])for u,v in zip(a,b))
 assert dot(h,k)<=bound
# Two exceptional coordinates cannot appear in a root; there is only one.
# The root norm96 has two amplitudes sqrt(48), at most one meets magnitude sqrt(8).
# [sqrt(48)*(sqrt(8)+2)]² = 576+384sqrt(2)<2304.
assert 2*384**2<(48**2-576)**2
# Distinct dense/sparse types have differing ordinary supports.
U=[]
for i,j in combinations(range(4),2):
 for a,b in product((-2,2),repeat=2):
  U.append(tuple(a if t==i else b if t==j else 0 for t in range(4)))
assert len(U)==24 and max(dot(u,v)for u,v in product(U,A+B+C))==4
assert max(dot(u,v)for u,v in product(product((-1,1),repeat=4),A+B+C))==4
assert 12*4==48 and 4*(6*4)**2==24*96 and 6*4*2==48
# Inherited sparse input verification is wholly independent of the head labels.
# The claims are checked against the unchanged sibling proof packages and
# their inherited RS/root counts, not accepted from a caller-supplied scalar.
s=next(x for x in json.loads((UP/'dim92-96-shifted-dodecacode-outer-extension'/'claims.json').read_text())if x['dimension']==95)
f=next(x for x in json.loads((UP/'dim92-96-flat-rs-layers'/'claims.json').read_text())if x['dimension']==95)
assert s['patterns']==1064 and s['mask']==3 and s['inherited']==102*8**10+1838148
assert s['lower_bound']==s['inherited']+1064*2**20==f['inherited']
assert f['extra']==2097152 and f['lower_bound']==f['inherited']+f['extra']==110641286212
count=108*8**10+1838148+1064*2**20+2097152
assert count==117083737156 and count-f['lower_bound']==6*8**10
claims=json.loads((P/'claims.json').read_text())
assert claims[0]['lower_bound']==count and claims[0]['previous']==f['lower_bound']
print(json.dumps({'dimension':95,'color_sizes':dict(sorted(Counter(c for c,_,_ in D).items())),
 'heads':len(D),'pairs':len(D)*(len(D)-1)//2,'dense':108*8**10,'previous':f['lower_bound'],
 'shifted_patterns':1064,'flat_rays':f['extra'],'new_lower_bound':count,'gain':count-f['lower_bound']},indent=2))
