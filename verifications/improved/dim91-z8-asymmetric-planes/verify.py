from itertools import product,combinations
from collections import Counter
from pathlib import Path
import hashlib,json,sys,subprocess
P=Path(__file__).resolve().parent
UP=P.parent if len(sys.argv)==1 else Path(sys.argv[1]).resolve()
p=P/'z8_generator.json'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
G=json.loads(p.read_text())['generator']
assert len(G)==12 and all(len(v)==24 and v[i]%2==1 and not any(v[:i])for i,v in enumerate(G))
assert len({(a*G[0][0]%8,(a*G[0][1]+b*G[1][1])%8)for a,b in product(range(8),repeat=2)})==64
A=[tuple(1 if ((q>>2)^(q&1)*x^((q>>1)&1)*y)&1==0 else -1 for x,y in product(range(2),repeat=2))for q in range(8)]
B=[(-a[0],*a[1:])for a in A]
C=[tuple((2 if q<4 else -2)if i==q%4 else 0 for i in range(4))for q in range(8)]
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
T=[A,B,C]
for t in T:
 assert len(set(t))==8 and all(dot(v,v)==4 for v in t)
for a,b in combinations(T,2):assert max(dot(x,y)for x,y in product(a,b))==2
D=[(p['color'],tuple(p['prefix']),tuple(p['head']))for p in json.loads((P/'heads.json').read_text())]
assert len(D)==len({(c,a)for c,a,_ in D})==24
SCALE=2
Q2=SCALE*SCALE
norm=lambda h:88*Q2+dot(h,h)
for c,a,h in D:
 assert 0<=c<=2 and len(a)==2 and all(0<=k<=7 for k in a)
 assert len(h)==3 and all(type(k)==int for k in h)and dot(h,h)<=8*Q2
 # Inherited parent homogeneous minimum24 means distinct same-prefix body dot<=48-8=40.
 assert 2*(40*Q2+dot(h,h))<=norm(h)
for (c,a,h),(e,b,k) in combinations(D,2):
 body=44 if c!=e else 48-sum(dot(T[c][u],T[c][v])for u,v in zip(a,b))
 raw=body*Q2+dot(h,k)
 assert raw<=0 or 4*raw**2<=norm(h)*norm(k)
U=[tuple(a if t==i else b if t==j else 0 for t in range(4))
   for i,j in combinations(range(4),2)for a,b in product((-2,2),repeat=2)]
assert len(U)==24 and max(dot(u,v)for u,v in product(U,A+B+C))==4
for w in (10,11):
 for c,a,h in D:
  assert 4*Q2*(4*w)**2<=8*w*norm(h)
# Replay the precise transferred 88/89D sparse layer, ignoring its dense count.
q=UP/'dim88-89-golay-variable-punctures'
s=json.loads(subprocess.run([sys.executable,str(q/'verify.py')],check=True,text=True,capture_output=True).stdout)
assert s['full_weight_11_patterns']==982 and s['old_weight_10_patterns']==11
sparse=982*1024*144+11*256*72; dense=24*8**10; total=dense+sparse
previous=22359805794
assert dense==25769803776 and sparse==145004544 and total==25914808320
claim=json.loads((P/'claims_base.json').read_text())[0]
assert claim['dimension']==91 and claim['previous']==previous and claim['lower_bound']==total
print(json.dumps({'dimension':91,'heads':len(D),'head_norm_squared':dict(sorted(Counter(dot(h,h)for c,a,h in D).items())),
 'color_counts':dict(sorted(Counter(c for c,_,_ in D).items())),
 'head_pairs':len(D)*(len(D)-1)//2,'dense':dense,'sparse':sparse,'previous':previous,'lower_bound':total,'gain':total-previous},indent=2))
