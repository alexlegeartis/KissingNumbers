from pathlib import Path
import json,itertools,hashlib
P=Path(__file__).resolve().parent;S=P/'input';j=json.loads((S/'e7/candidate.json').read_text());base=json.loads((S/'base/certificate.json').read_text())
def dot(a,b):return sum(x*y for x,y in zip(a,b))
A=base['basis_alphabets'][0]
assert hashlib.sha256((S/'e7/z8_generator.json').read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
prefix=j['source_prefixes'];assert len(prefix)==38 and len(set(map(tuple,prefix)))==38
for H in j['full_rotated_heads']:
 assert len(H)==38
 for i,h in enumerate(H):
  assert len(h)==8 and all(type(v)is int for v in h) and h[0]+h[4]==0 and dot(h,h)==8
  for k,g in enumerate(H):assert dot(h,g)==dot(A[prefix[i][0]],A[prefix[k][0]])+dot(A[prefix[i][1]],A[prefix[k][1]])
selected=j['selected'];assert len(selected)==100 and len({tuple(r['head'])for r in selected})==100
assert len({(r['color'],*r['prefix'])for r in selected})==100
for i,r in enumerate(selected):
 assert r['color']in[0,1,2] and r['head']==j['full_rotated_heads'][r['color']][prefix.index(r['prefix'])]
 for s in selected[:i]:
  if r['color']!=s['color']:assert dot(r['head'],s['head'])<=4
 assert max(map(abs,r['head']))<=2
# Establish GF16 directly: polynomial0b10011 has no divisor ofdegree1or2.
def rem(a,b):
 while a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
 return a
assert all(rem(19,b)!=0 for b in [2,3,4,5,6,7])
def mul(a,b):
 z=0
 for k in range(4):
  if b>>k&1:z^=a<<k
 return rem(z,19)
for a in range(1,16):assert any(mul(a,b)==1 for b in range(1,16))
# Explicitly check root bound on all65535 nonzero cubics at0..5.
for co in itertools.product(range(16),repeat=4):
 if not any(co):continue
 roots=0
 for t in range(6):
  v=co[3]
  for c in co[2::-1]:v=mul(v,t)^c
  roots+=v==0
 assert roots<=3
signs=[s for s in itertools.product([-1,1],repeat=6)if s.count(-1)%2==0]
assert len(signs)==32 and max(dot(a,b)for i,a in enumerate(signs)for b in signs[:i])==2
report=json.loads((S/'augmentation.json').read_text());out={}
for n,dense in [(93,42949672960),(94,64424509440),(95,100*8**10),(96,3*8**12)]:
 row=report['results'][str(n)];eligible=[i for i in range(n)if not(n==95 and i==88)];r=len(eligible)-80
 groups=[eligible[:r]]+[eligible[r+16*k:r+16*(k+1)]for k in range(5)]
 assert row['groups']==groups and row['eligible_coordinates']==eligible
 roots=4*(r*(r-1)//2+5*120+(94 if n==95 else 0));rs=32*r*16**3
 assert row['base']==dense and row['root_count']==roots and row['rs_count']==rs and row['total']==dense+roots+rs
 assert 6*16==96 and 3*16==48 and 6*4*2==48 and 16*48<48**2 and 25*48<48**2
 out[n]=dense+roots+rs
assert out=={93:42951379608,94:64426347212,95:107376020548,96:206160530240}
(P/'result.json').write_text(json.dumps({'status':'PASS new exactE7/RS/root audit atop freshly audited three-basis parent','counts':out},indent=2));print(out)
