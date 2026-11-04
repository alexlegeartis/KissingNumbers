"""Independent standard-library replay of factor coverage and finite obstructions.
Uses no NumPy, SymPy, or search-module imports. Raises on invalid evidence.
"""
from pathlib import Path
from collections import Counter
import json,hashlib,time,re
OUT=Path(__file__).resolve().parent

def require(ok,message):
 if not ok:raise ValueError(message)
def rem(a,b):
 while a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
 return a
def gcd(a,b):
 while b:a,b=b,rem(a,b)
 return a
def square(a):return sum(((a>>i)&1)<<(2*i) for i in range(a.bit_length()))
def irreducible(p):
 d=p.bit_length()-1;x=rem(2,p);v=x
 for i in range(1,d+1):
  v=rem(square(v),p)
  if i<=d//2 and gcd(v^x,p)!=1:return False
 return v==x
def polyproduct(a,b,q):
 v=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%q
 return v
def packed(p):return sum(x<<i for i,x in enumerate(p))
def weight(w,q):return sum(0 if x==0 else q//2 if x==q//2 else q//4 for x in w)
def word(row,c):
 q=row['q'];n=row['n'];g=row['generator_polynomial'];k=row['rank']
 require(len(c)==k and all(type(x)is int and 0<=x<q for x in c),'bad message')
 v=[0]*n
 for i,a in enumerate(c):
  for j,b in enumerate(g):v[i+j]=(v[i+j]+a*b)%q
 return v
def witness(row,w):
 v=word(row,w['message']);require(v==w['word'],'word mismatch')
 require(any(v),'zero witness');require(weight(v,row['q'])==w['base_weight'],'weight mismatch')
 return v

def inconsistent(constraints,k):
 # Independent elementary equation elimination over Z4; unit pivots first.
 A=[c['message'][:]+[c['required_extension_value']]for c in constraints];rank=0
 for j in range(k):
  pivot=next((i for i in range(rank,len(A))if A[i][j]%2),None)
  if pivot is None:continue
  A[rank],A[pivot]=A[pivot],A[rank];inverse=A[rank][j]
  A[rank]=[v*inverse%4 for v in A[rank]]
  for i in range(len(A)):
   if i==rank:continue
   factor=A[i][j];A[i]=[(a-factor*b)%4 for a,b in zip(A[i],A[rank])]
  rank+=1
 B=[]
 for r in A[rank:]:
  require(all(x%2==0 for x in r[:-1]),'non-even residual')
  if r[-1]%2:return True
  B.append(sum((x//2)<<j for j,x in enumerate(r)))
 pivots={}
 for r in B:
  for j in range(k):
   if not(r>>j&1):continue
   if j in pivots:r^=pivots[j]
   else:pivots[j]=r;break
  else:
   if r==1<<k:return True
 return False

def main():
 start=time.monotonic();data=json.loads((OUT/'composite_screen.json').read_text())
 ROOT=OUT.parents[2]
 for filename,digest in data['baseline_source_hashes'].items():require(hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==digest,'baseline source changed')
 references={int(row[0]):int(row[1].replace(' ',''))for row in json.loads((ROOT/'research/group_search_20261104/baselines.json').read_text())['rows']}
 for c in json.loads((ROOT/'discovery/benchmarks/catalog.json').read_text())['comparisons']:references[c['dimension']]=max(references.get(c['dimension'],0),c['cardinality'])
 for line in (ROOT/'README.md').read_text().splitlines():
  m=re.match(r'\| (\d+) \|.*?\| \*\*([0-9 \u202f]+)\*\* \|',line)
  if m:references[int(m[1])]=max(references.get(int(m[1]),0),int(re.sub(r'\s','',m[2])))
 allowed={(n,q)for q in(4,8,16)for n in range(3,48,2)if 25<=(n+1)*q//2<=96}
 actual=[(f['n'],f['q'])for f in data['families']]
 require(set(actual)==allowed and len(actual)==len(allowed),'family coverage mismatch')
 rows={(r['n'],r['q'],r['factor_subset']):r for r in data['examined']}
 require(len(rows)==len(data['examined']),'duplicate module')
 expected=set();pruned=0;total=0
 for family in data['families']:
  n=family['n'];q=family['q'];factors=family['binary_factors'];p=[1]
  reference_dimension=max(d for d in references if d<=(n+1)*q//2)
  require(family['baseline_dimension']==reference_dimension and family['baseline']==references[reference_dimension],'baseline snapshot mismatch')
  require(25<=(n+1)*q//2<=96 and q in(4,8,16)and n%2==1,'family range')
  for f in factors:
   require(irreducible(packed(f)),'nonirreducible factor');p=polyproduct(p,f,2)
  require(p==[1]+[0]*(n-1)+[1],'incomplete factorization')
  for subset in range(1<<len(factors)):
   g=[1]
   for i,f in enumerate(factors):
    if subset>>i&1:g=polyproduct(g,f,2)
   k=n-len(g)+1
   if k<=0:continue
   total+=1
   if q**k*2**(q//2)<=family['baseline']:pruned+=1;continue
   key=(n,q,subset);expected.add(key);r=rows[key]
   require(r['rank']==k and r['cardinality']==q**k and r['binary_length']==(n+1)*q//2,'module parameters')
   require(r['target_distance']==(r['binary_length']+3)//4,'wrong target')
   require([x%2 for x in r['generator_polynomial']]==g,'wrong divisor')
   require(r['generator_polynomial'][0]%2==1 and r['generator_polynomial'][-1]==1,'nonfree generator')
   require(polyproduct(r['generator_polynomial'],r['complementary_polynomial'],q)==[q-1]+[0]*(n-1)+[1],'lift identity')
 require(expected==set(rows),'module coverage mismatch');require(pruned==data['capacity_pruned'],'pruning mismatch')
 classifications={}
 for key,r in rows.items():
  q=r['q'];target=r['target_distance'];v=witness(r,r['best'])
  for scalar,w in r['scalar_witnesses'].items():
   v=witness(r,w);a=int(scalar)*sum(v)%q
   require(a==w['extension_value']and weight(v,q)+weight([a],q)==w['extended_weight']<target,'scalar rejection invalid')
  if r['status']=='all_extra_coordinate_values_excluded':
   require(r['best']['base_weight']+q//2<target,'universal-coordinate obstruction invalid');classifications[key]='every_extra_coordinate_value'
  elif r['status']=='all_scalar_sum_extensions_excluded':
   require(set(r['scalar_witnesses'])==set(map(str,range(q))),'missing scalar');classifications[key]='every_scalar_sum_extension'
  else:classifications[key]='unresolved'
 refine=json.loads((OUT/'refine_residues.json').read_text())
 require(refine['input_sha256']==hashlib.sha256((OUT/'composite_screen.json').read_bytes()).hexdigest(),'stale refinement input')
 for w in refine['results']:
  key=(w['n'],w['q'],w['factor_subset']);r=rows[key];v=witness(r,w);target=r['target_distance'];q=r['q']
  if w['excludes_every_extra_coordinate_value']:
   require(weight(v,q)+q//2<target,'refinement obstruction invalid');classifications[key]='every_extra_coordinate_value'
  else:
   for s in w['scalar_witnesses']:require(weight(v,q)+weight([int(s)*sum(v)%q],q)<target,'refinement scalar obstruction invalid')
   if set(w['scalar_witnesses'])|set(r['scalar_witnesses'])==set(map(str,range(q))):classifications[key]='every_scalar_sum_extension'
 linear=json.loads((OUT/'linear31_obstruction.json').read_text())
 for proof in linear['results']:
  if not proof['inconsistent']:continue
  key=(proof['n'],proof['q'],proof['factor_subset']);r=rows[key]
  require(proof['q']==r['q']==4 and bool(proof['constraints']),'linear proof domain')
  for c in proof['constraints']:
   v=witness(r,c);require(weight(v,4)==14 and r['target_distance']==16 and c['required_extension_value']==2,'linear necessity invalid')
  require(inconsistent(proof['constraints'],r['rank']),'linear equations consistent')
  multipliers=proof['contradiction_multipliers'];require(multipliers is not None,'missing elementary contradiction')
  require(all(sum(a*c['message'][j]for a,c in zip(multipliers,proof['constraints']))%4==0 for j in range(r['rank'])),'left side not zero')
  require(sum(a*c['required_extension_value']for a,c in zip(multipliers,proof['constraints']))%4==2,'right side not 2')
  classifications[key]='every_linear_extension'
 result=dict(status='PASS',modules_total=total,modules_capacity_pruned=pruned,modules_examined=len(rows),
  classification=dict(Counter(classifications.values())),unresolved=[list(k)for k,v in classifications.items()if v=='unresolved'],
  elapsed_seconds=time.monotonic()-start,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (OUT/'independent_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
