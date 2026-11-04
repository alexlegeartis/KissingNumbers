"""Bounded exact screen of all coprime cyclic-factor products; no discovery claims.
Only these local files are written. Run from repository root.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import hashlib, json, re, sys, time, warnings
import numpy as np
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from discovery.engines.rings import _solve2
warnings.filterwarnings('ignore')
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent

def lift(g,h,n,q):
 g=np.array(g,dtype=np.int64);h=np.array(h,dtype=np.int64)
 a,b=len(g)-1,len(h)-1
 target=np.zeros(n+1,dtype=np.int64);target[0]=-1;target[-1]=1
 current=2
 while current<q:
  residual=target-np.convolve(g,h)
  if np.any(residual%current): raise ValueError('bad factorization')
  M=np.zeros((n,n),dtype=np.int64)
  for i in range(a):M[i:i+b+1,i]=h
  for i in range(b):M[i:i+a+1,a+i]=g
  correction=_solve2(M,(residual//current)[:-1])
  g[:-1]+=current*correction[:a];h[:-1]+=current*correction[a:]
  current*=2;g%=current;h%=current
  if not np.array_equal(np.convolve(g,h)%current,target%current):raise ValueError('lift identity')
 return g.tolist(),h.tolist()

def convolution(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]^=x*y
 return c

def word(G,coef,q):return (np.array(coef,dtype=np.int64)@G%q).tolist()
def weight(w,q):return sum(0 if x==0 else q//2 if x==q//2 else q//4 for x in w)

def screen(G,q,target):
 k,n=G.shape
 best=None; witnesses={}; trials=0
 # Exact two-support messages, then bounded three-support messages.
 for support_size,values in [(1,range(1,q)),(2,range(1,q)),(3,sorted({1,2,q//2}))]:
  for support in combinations(range(k),support_size):
   for coeff in product(values,repeat=support_size):
    c=[0]*k
    for i,v in zip(support,coeff):c[i]=v
    w=word(G,c,q);wt=weight(w,q);trials+=1
    if best is None or wt<best['base_weight']:best=dict(message=c,word=w,base_weight=wt)
    for s in range(q):
     if str(s) in witnesses:continue
     extended=(s*sum(w))%q
     if wt+weight([extended],q)<target:
      witnesses[str(s)]=dict(message=c,word=w,base_weight=wt,extension_value=extended,extended_weight=wt+weight([extended],q))
    if wt+q//2<target:
     return dict(status='all_extra_coordinate_values_excluded',best=best,scalar_witnesses=witnesses,messages_checked=trials)
    if len(witnesses)==q:
     return dict(status='all_scalar_sum_extensions_excluded',best=best,scalar_witnesses=witnesses,messages_checked=trials)
 return dict(status='unresolved_after_bounded_message_scan',best=best,scalar_witnesses=witnesses,messages_checked=trials)

def main():
 start=time.monotonic();x=sp.Symbol('x')
 baselines={int(row[0]):int(row[1].replace(' ','')) for row in json.loads((ROOT/'research/group_search_20261104/baselines.json').read_text())['rows']}
 for c in json.loads((ROOT/'discovery/benchmarks/catalog.json').read_text())['comparisons']:
  baselines[c['dimension']]=max(baselines.get(c['dimension'],0),c['cardinality'])
 for line in (ROOT/'README.md').read_text().splitlines():
  match=re.match(r'\| (\d+) \|.*?\| \*\*([0-9 \u202f]+)\*\* \|',line)
  if match:
   dimension=int(match[1]); cardinality=int(re.sub(r'\s','',match[2]))
   baselines[dimension]=max(baselines.get(dimension,0),cardinality)
 rows=[];families=[];pruned=0
 for q in (4,8,16):
  for n in range(3,48,2):
   N=(n+1)*q//2
   if not 25<=N<=96:continue
   f=sp.Poly(x**n-1,x,modulus=2)
   factors=[[int(p.nth(i))%2 for i in range(p.degree()+1)]for p,m in sp.factor_list(f)[1]]
   total=0
   for subset in range(1<<len(factors)):
    degree=sum(len(v)-1 for i,v in enumerate(factors) if subset>>i&1)
    rank=n-degree
    if rank<=0:continue
    total+=1
    cardinality=q**rank
    # At most 2^(local binary length) distinct repeated binary masks exist.
    max_repeated_cosets=2**(q//2)
    baseline_dimension=max(d for d in baselines if d<=N)
    baseline=baselines[baseline_dimension]
    if cardinality*max_repeated_cosets<=baseline:
     pruned+=1;continue
    g=h=[1]
    for i,v in enumerate(factors):
     if subset>>i&1:g=convolution(g,v)
     else:h=convolution(h,v)
    g,h=lift(g,h,n,q)
    G=np.zeros((rank,n),dtype=np.int64)
    for i in range(rank):G[i,i:i+len(g)]=g
    result=screen(G,q,(N+3)//4)
    row=dict(n=n,q=q,binary_length=N,rank=rank,cardinality=cardinality,baseline=baseline,baseline_dimension=baseline_dimension,
       max_repeated_cosets=max_repeated_cosets,factor_subset=subset,binary_factors=factors,
       generator_polynomial=g,complementary_polynomial=h,target_distance=(N+3)//4,**result)
    rows.append(row)
   families.append(dict(n=n,q=q,binary_length=N,all_nonzero_modules=total,binary_factor_degrees=[len(v)-1 for v in factors],binary_factors=factors,baseline_dimension=max(d for d in baselines if d<=N),baseline=baselines[max(d for d in baselines if d<=N)]))
 result=dict(scope='Every monic divisor of x^n-1, odd n, q in 4/8/16, 25<=binary length<=96; capacity pruning pertains only to unions using all repeated binary masks, at most 2^(q/2) translates. It does not prune general nonconstant masks or other assembly.',
  status='bounded exact rejection search; unresolved is not certified',families=families,capacity_pruned=pruned,examined=rows,
  summary=dict(Counter(r['status'] for r in rows)),elapsed_seconds=time.monotonic()-start,
  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),baseline_source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'README.md',ROOT/'discovery/benchmarks/catalog.json',ROOT/'research/group_search_20261104/baselines.json')})
 (OUT/'composite_screen.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k]for k in ('summary','capacity_pruned','elapsed_seconds')}))
 print('families',len(families),'examined',len(rows))
 for r in rows:
  if r['status']=='unresolved_after_bounded_message_scan':print('UNRESOLVED',r['n'],r['q'],r['rank'],r['best']['base_weight'],sorted(r['scalar_witnesses']))
if __name__=='__main__':main()
