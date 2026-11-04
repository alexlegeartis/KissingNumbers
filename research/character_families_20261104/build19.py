"""Rational variant of Li's current19D12460 control, from the printed formulas.

The diagonal tail factor817/1000 replaces sqrt(2/3), and the mark is747/1000.
All points are independently
checked through the existing unequal-norm backend; this is not a new record.
"""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
OUT=Path(__file__).resolve().parent
GEN=[(1,8,9,12,16,17,18,19),(2,10,11,14,15,17,18,20),(3,7,9,13,15,16,17,19),(4,7,8,10,12,15,16,19),(5,10,12,13,15,16,17,18),(6,7,8,9,10,13,16,18),(1,4,7,9),(1,5,6,18,20),(1,3,12,15,20),(1,10,13,19,20),(1,3,5,6,7,13,14,15,18,20),(2,4,6,7,8,13,14,16,17,18)]
REPS=[1134,2453,4337,5590,6500,7434,8363,9413,10359,11600,12378,15592,17037,18346,19224,20342,21301,24455,26184,27571,29399,30393,31243,32556,66280,69210,70169,71287,72389,73698,74307,75620,76758,77752,79580,80679,82021,82955,84153,85406,87290,88321,90503,93493,94271,95512,96682,97732,986,1972,2822,3617,4907,5644,6846,7888,9073,10015,10980,12227,13184,13991,15196,16178,16414,17520,18626,19941,20902,21633,22579,23645,25084,26002,26729,27982,28740,30051,30872,31990,65842,66581,67751,68809,70083,71085,71967,72760,73936,75255,75788,76898,77857,78927,80308,81043,82678,83921,84835,85773,86862,87840,88978,89781,90717,92026,92801,93935,95205,96139,96880,98135]
T=[0,3,6,8];H=[j for j in range(20) if j not in T];Tm=sum(1<<j for j in T)
R=np.array([[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]],dtype=np.int64)
def span(g):
 out=[0]
 for b in g:out += [x^b for x in out]
 return sorted(out)
def mask(s):return sum(1<<(j-1) for j in s)
def codes():
 E=span(list(map(mask,GEN)))
 # Nullspace by exact binary elimination, then enumerate only2^8 words.
 piv={}
 for row in map(mask,GEN):
  v=row
  while v:
   i=v.bit_length()-1
   if i in piv:v^=piv[i]
   else:piv[i]=v;break
 free=[j for j in range(20) if j not in piv];dual_basis=[]
 for j in free:
  v=1<<j
  for i in sorted(piv):
   if (v&piv[i]).bit_count()%2:v|=1<<i
  dual_basis.append(v)
 dual=span(dual_basis)
 if any((x&g).bit_count()%2 for x in dual for g in map(mask,GEN)):raise ValueError('dual')
 octads=[w for w in dual if w.bit_count()==8]
 if (len(E),len(dual),len(octads))!=(4096,256,130):raise ValueError('code parameters')
 tetrads=[w for w in E if w.bit_count()==4]
 U0=span([w^Tm for w in tetrads if w!=Tm])
 if len(U0)!=16:raise ValueError('tetrad group')
 return E,octads,U0

def build(mark=747):
 E,octads,U0=codes();body=[]
 def append(row):
  if sum(row[j] for j in T):return
  tail=R@np.array([row[j] for j in T]);assert np.all(tail%2==0)
  body.append([row[j] for j in H]+list(map(int,tail//2)))
 for i,j in combinations(range(20),2):
  for a,b in product((-2,2),repeat=2):
   row=[0]*20;row[i]=a;row[j]=b;append(row)
 for o in octads:
  positions=[j for j in range(20) if o>>j&1]
  for signs in product((-1,1),repeat=8):
   if sum(s<0 for s in signs)%2!=1:continue
   row=[0]*20
   for j,s in zip(positions,signs):row[j]=s
   append(row)
 W=sorted({c^u for c in REPS for u in U0});Wset=set(W)
 if len(W)!=1792 or any(c not in E for c in W):raise ValueError('W')
 diffs=[d for d in E if (d&Tm).bit_count()==1 and (d&~Tm).bit_count()==4]
 additions=[];metadata=[]
 for c in W:
  q=np.array([1-2*(c>>j&1) for j in range(20)]);t=R@q[T]//2
  row=[1000*int(q[j]) for j in H]
  if (c&Tm).bit_count()==2:
   forbidden=0
   for d in diffs:
    if c^d in Wset:forbidden|=d
   marks=[j for j in H if not forbidden>>j&1]
   if len(marks)!=1:raise ValueError('free coordinate')
   row[H.index(marks[0])]=mark*int(q[marks[0]]);row+=list(map(int,1000*t))
  else:row+=list(map(int,817*t))
  additions.append(row);metadata.append(c)
 if len(body)!=10668:raise ValueError('body count')
 # Sign action of U0 is scalar on T, so it descends to all three tail axes.
 g=[]
 tetrads=[w for w in E if w.bit_count()==4 and w!=Tm]
 for w in tetrads:
  u=w^Tm;g.append([u>>j&1 for j in H]+[1]*3)
 return body,additions,metadata,g

def main():
 start=time.monotonic();body,additions,words,g=build();rows=body+additions
 artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':rows,'radical_coordinates':[[0]*19 for _ in rows]})
 context=RunContext(ROOT,{'max_pairs':80000000,'max_materialized_points':15000})
 report=verify_artifact(artifact,context);artifact.verification=report;artifact.facts=report.details['facts']
 (OUT/'control19_artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 program={'schema_version':1,'name':'Known19D12460 rational variant with strict auxiliary gap',
          'description':'Printed Li coset choices, rational tail817/1000 and mark747/1000. Known record count, no held-out discovery.',
          'nodes':[{'id':'sphere','op':'geometry.quadratic_code','params':artifact.data}], 'outputs':['sphere']}
 (OUT/'control19_program.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
 # Every pair involving an addition has a strictly positive angular margin.
 # All raw norms take only three values. Exact maxima per norm class reduce
 # the final large-integer gap decision to five comparisons.
 values=np.array(rows,np.int64);norms=np.sum(values*values,axis=1)
 maxima={}
 for at in range(len(body),len(rows),64):
  dots=values[at:at+64]@values.T
  for offset,products in enumerate(dots):
   i=at+offset
   for other in (8,18002467,19558009):
    ids=np.flatnonzero(norms==other);ids=ids[ids!=i]
    if len(ids):
     k=tuple(sorted((int(norms[i]),other)));maxima[k]=max(maxima.get(k,-10**20),int(products[ids].max()))
 gap=[]
 for (a,b),dot in sorted(maxima.items()):
  if dot>0 and (30000*dot)**2>14999**2*a*b:raise ValueError('Strict auxiliary gap failed')
  gap.append({'norms':[a,b],'maximum_raw_dot':dot,'cosine_upper_bound':[14999,30000]})
 (OUT/'control19_source.json').write_text(json.dumps({'source':'https://hexagonmath.org/2610.00248v1','source_polynomial19':'Ho twelve generators printed in Section4','source_W':'Table5:112 representatives','new_choice':'diagonal tail817/1000 in place of sqrt(2/3), axis mark747/1000 in place of373/500','body_count':len(body),'addition_count':len(additions),'generators':g,'addition_words':words,'facts':artifact.facts,'verification':report.to_dict(),'auxiliary_gap':gap,'budget':context.budget,'seconds':time.monotonic()-start,'new_sota':False,'scope':'Known lower-bound count12460, rational-coordinate variant; supplied printed W table, no held-out discovery'},indent=2)+'\n')
 print(json.dumps({'facts':artifact.facts,'pairs_checked':report.details['pairs_checked'],'seconds':time.monotonic()-start}),flush=True)
if __name__=='__main__':main()
