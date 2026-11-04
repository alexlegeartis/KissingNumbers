"""Formula reconstruction of the known5730-point17D layered control.

The established family ceiling is unchanged. This is not a record search.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import json,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.geometry import integer_code

def main():
 out=Path(__file__).resolve().parent;start=time.monotonic()
 rm1=set()
 for f in range(16):
  for c in (0,1):rm1.add(sum((((x&f).bit_count()%2)^c)<<x for x in range(16)))
 if len(rm1)!=32:raise ValueError('RM1 dimension')
 rm2=[b for b in range(65536) if all((b&r).bit_count()%2==0 for r in rm1)]
 m0=sum(((((x>>0)&1)*((x>>1)&1))^(((x>>2)&1)*((x>>3)&1)))<<x for x in range(16))
 supports=sorted(s for r in rm1 if (s:=m0^r).bit_count()==6)
 if len(supports)!=16:raise ValueError('Bent quadratic hexad orbit')
 kl=sorted(b for b in rm2 if (b&m0).bit_count()%2==0);w4=[b for b in kl if b.bit_count()==4]
 if (len(rm2),len(kl),len(w4))!=(2048,1024,60):raise ValueError('TierB code parameters')
 fibers={};reps=[]
 for b in kl:
  if b not in fibers:
   for r in rm1:fibers[b^r]=len(reps)
   reps.append(b)
 adjacent=set()
 for b in kl:
  for w in w4:
   a,c=fibers[b],fibers[b^w]
   if a!=c:adjacent.add((min(a,c),max(a,c)))
 six=next(s for s in combinations(range(len(reps)),6) if all((a,b) not in adjacent for a,b in combinations(s,2)))
 plus=sorted(b for b in kl if fibers[b] in six);minus=sorted(b^w4[0] for b in plus)
 if len(plus)!=192 or set(plus)&set(minus):raise ValueError('Disjoint192-word tiers')
 rows=[]
 # Squared norm72 in metric diag(1,...,1,2).
 for i,j in combinations(range(16),2):
  for a,b in product((-6,6),repeat=2):
   row=[0]*17;row[i]=a;row[j]=b;rows.append(row)
 for support in sorted(r for r in rm1 if r.bit_count()==8):
  positions=[j for j in range(16) if support>>j&1]
  for signs in range(256):
   if signs.bit_count()%2!=1:continue
   row=[0]*17
   for k,j in enumerate(positions):row[j]=-3 if signs>>k&1 else 3
   rows.append(row)
 for support in supports:
  positions=[j for j in range(16) if support>>j&1]
  for signs in range(64):
   if signs.bit_count()%2!=1:continue
   for height in (-3,3):
    row=[0]*16+[height]
    for k,j in enumerate(positions):row[j]=-3 if signs>>k&1 else 3
    rows.append(row)
 for words,height in ((plus,2),(minus,-2)):
  for word in words:rows.append([-2 if word>>j&1 else 2 for j in range(16)]+[height])
 for height in (-6,6):rows.append([0]*16+[height])
 context=RunContext(ROOT,{'max_pairs':20000000})
 params={'coordinates':rows,'weights':[1]*16+[2],'norm_squared':72}
 artifact=integer_code([],params,context)
 if artifact.facts['cardinality']!=5730 or artifact.facts['dimension']!=17:raise ValueError('Control count/rank')
 program={'schema_version':1,'name':'Formula-reconstructed known17D5730 control','description':'Known layered construction; supplied runtime coordinate artifact has formula-generation lineage, not held-out discovery.','nodes':[{'id':'sphere','op':'geometry.integer_code','params':params}],'outputs':['sphere']}
 (out/'control17_program.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
 (out/'control17_artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 report={'facts':artifact.facts,'verification':artifact.verification.to_dict(),'budget':context.budget,
         'heights':dict(Counter(row[-1] for row in rows)),'quadratic_support_seed':m0,'independent_fibers':list(six),
         'seconds':time.monotonic()-start,'new_sota':False,'scope':'Known5730 reconstruction in a previously closed family. Not a new record or a new search of that family.'}
 (out/'control17_summary.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'facts':artifact.facts,'heights':report['heights'],'seconds':report['seconds']}),flush=True)
if __name__=='__main__':main()
