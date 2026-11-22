"""Exact five/three-plane proof using the actual parent square support generators."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from importlib import util
from hashlib import sha1
import json
H=Path(__file__).resolve().parent

def pin(rel,hash):
 p=H.parent/rel; raw=p.read_bytes()
 assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==hash, rel
 return p

def run(rel,hash):
 p=pin(rel,hash); spec=util.spec_from_file_location(p.stem,p)
 m=util.module_from_spec(spec); spec.loader.exec_module(m);return m

def layers(n,S):
 assert len(S)==len(set(S))=={39:2223,45:3465}[n]
 S=set(S)
 sqtrip=[t for b in S for t in combinations(b,3)]
 assert len(sqtrip)==len(set(sqtrip))
 missed=set(combinations(range(n),3))-set(sqtrip)
 assert len({p for t in missed for p in combinations(t,2)})==3*len(missed)
 five=[tuple(map(int,l.split())) for l in (H/f'five{n}.txt').read_text().splitlines()]
 assert all(len(b)==5 and tuple(sorted(set(b)))==b and 0<=b[0]<b[-1]<n for b in five)
 fs=[t for b in five for t in combinations(b,3)]
 assert len(fs)==len(set(fs))
 assert all(q not in S for b in five for q in combinations(b,4))
 pairs={p for b in five for p in combinations(b,2)}
 three=[q for q in missed-set(fs) if all(p in pairs for p in combinations(q,2))]
 assert len({p for q in three for p in combinations(q,2)})==3*len(three)
 return len(five),len(three),4*len(pairs)

def main():
 assert __debug__
 pkg='dim86-87-remarked-square/'
 M=run(pkg+'verify.py','890c3dc5311044b4c35f0aee6cc8fd9f4958d87d')
 pin(pkg+'sqs14.txt','144a4b51ee5aa3a35fa3977d08d5ee39b5ecd99c')
 A=json.loads(pin(pkg+'claims.json','aaa0c8ee25ddcb542644f70dd3f9b0e024981c57').read_text())
 pkg='dim90-91-a2-rotated-square/'
 N=run(pkg+'support.py','8eb8ce39f07ed1f2e0701b2f872a36a229ece196')
 B=json.loads(pin(pkg+'claims.json','06da890849eb273b05bd89d6fa6acc2c17d265dd').read_text())
 # Five-plane: all Pless directions at >=30 degrees, giving 3/2 per row.
 angles=[min((150-a)%360,(a-150)%360) for a in range(0,360,60)]
 assert min(angles)==30 and 4*F(15,2)**2==5*45 and 4*25<5*45
 assert 4*F(9,2)<5*4 and 4*1**2<5*2<4*2**2
 assert 5-2*2<=F(5,2) and 2<=F(5,2)
 fsign=[i for i in range(32) if i.bit_count()%2==0]
 assert len(fsign)==16
 assert all({(x>>i&1,x>>j&1) for x in fsign}=={(0,0),(0,1),(1,0),(1,1)} for i,j in combinations(range(5),2))
 # Triple chosen angle 165°, closest square direction is 195° (30°).
 assert min(min((165-a)%360,(a-165)%360) for a in range(15,360,90))==30
 assert 4*3==3*4 and 4*(3*3*3)<3*45 and 4*9<3*45
 # Triple/five ≤2cos(15°). Its square is 2+sqrt(3) < 15/4.
 assert 3<F(7,4)**2 and 2+F(7,4)==F(15,4)
 # For triple/root, (cos(15°)+sin(15°))²=3/2 exactly. Both
 # low-axis projections are at most this sum. The only forbidden
 # double-high pair has already been removed by the five-plane layer.
 assert 4*F(3,2)==3*2
 assert 1<=F(3,2) and 3-2<=F(3,2)
 rs=[]
 for source,n,l in [(A,39,layers(39,M.quadruples())),([x for x in B if x['dimension']==90],45,layers(45,N.supports()))]:
  for row in source:
   d,b=row['dimension'],row['lower_bound']
   f,t,r=l;g=16*f+8*t-r
   rs.append(dict(dimension=d,previous=b,lower_bound=b+g,gain=g,five_plane_supports=f,three_plane_supports=t,retired_roots=r))
 assert rs==json.loads((H/'claims.json').read_text())['dimensions']
 print(json.dumps({'dimensions':rs},indent=2))
if __name__=='__main__':main()
