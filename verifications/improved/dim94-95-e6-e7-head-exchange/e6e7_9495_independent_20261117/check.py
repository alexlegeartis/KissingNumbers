from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
A=[tuple((-1)**((a//4)^((a%2)*x)^(((a//2)%2)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])for a in range(8)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
results=[]
for n,name,count,extra in [(94,'e7_to_e6_94_20261117',62,1837772),(95,'e7_mixed_sections95_20261117',102,1838148)]:
 S=P/name;j=json.loads((S/'candidate.json').read_text());assert hashlib.sha256((S/'z8_generator.json').read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
 rows=j['selected'];assert len(rows)==count and len({tuple(r['head'])for r in rows})==count and len({(r['color'],*r['prefix'])for r in rows})==count
 for i,r in enumerate(rows):
  h=r['head'];assert r['color']in range(3)and len(r['prefix'])==2 and all(a in range(8)for a in r['prefix'])
  assert len(h)==8 and all(type(x)is int for x in h)and dot(h,h)==8 and h[0]+h[4]==0 and max(map(abs,h))<=2
  if n==94:assert h[1]==-h[0] and abs(h[0])<=1
  for s in rows[:i]:
   if r['color']==s['color']:
    a,b=r['prefix'];c,d=s['prefix'];assert dot(h,s['head'])==dot(A[a],A[c])+dot(A[b],A[d])
   else:assert dot(h,s['head'])<=4
 # Coordinate transform squared norm and all used coordinate bounds are exact.
 for r in rows:
  h=r['head'];assert (3*h[0]**2+sum(h[k]**2 for k in [2,3,5,6,7])if n==94 else 2*h[0]**2+sum(h[k]**2 for k in [1,2,3,5,6,7]))==8
 rs=32*14*16**3;roots=2*14*13+5*2*16*15+(4*94 if n==95 else 0);assert extra==rs+roots
 total=count*8**10+extra;assert total=={94:66573830860,95:109523504196}[n]
 results.append({'dimension':n,'head_count':count,'dense_count':count*8**10,'RS_root_count':extra,'total':total,'candidate_sha256':hashlib.sha256((S/'candidate.json').read_bytes()).hexdigest(),'status':'PASS new exact head geometry; source and RS/root proofs inherited from completed audit'})
(P/'result.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
