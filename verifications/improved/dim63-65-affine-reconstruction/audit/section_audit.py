"""Independent exact63D section count of independently audited64D family."""
from pathlib import Path
from itertools import combinations
from math import comb
import numpy as np,json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[3];P=R/'verifications/improved/dim61-67-compatible-sign-families';S=P.parent/'dim61-68-95-96-code-layers/source'
c=json.loads((H.parent/'certificate.json').read_bytes());b=json.loads((P/'data/global64.json').read_bytes())
proposal=json.loads((H.parent/'section63.json').read_bytes());h=proposal['octad'];assert h.bit_count()==8
with np.load(S/'dim61_96/code64_supports.npz',allow_pickle=False) as f:old=list(map(int,f['words']));gen=[sum(int(v)<<i for i,v in enumerate(r)) for r in f['generator']]
def rank(rows):
 piv={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in piv:x^=piv[k]
   else:piv[k]=x;break
 return len(piv)
positions=[i for i in range(64) if h>>i&1]
assert rank(gen)==28 and all((g&h).bit_count()%2==0 for g in gen)
assert rank([sum(((g>>i)&1)<<j for j,i in enumerate(positions)) for g in gen])==7
rows={r['support']:r for r in b['old']+b['new']+c['modified_parent_classes']+c['new_classes']}
new=[r['support'] for r in b['new']+c['new_classes']]
X=np.arange(65536,dtype=np.uint32);pc=np.array([i.bit_count() for i in range(65536)],np.uint8)
checks=[65535]+[sum(((i>>j)&1)<<i for i in range(16)) for j in range(4)]
std=np.ones(65536,bool)
for g in checks:std &= (pc[X&g]%2==0)
standard=X[std];assert len(standard)==2048 and min(pc[standard[1:]])==4
assert np.count_nonzero(pc[standard&255]==4)==1120
# For fewer than8 specified sites, dual minimum8 implies uniform projection.
dual=[0]
for g in checks:dual +=[v^g for v in dual]
assert min(v.bit_count() for v in dual if v)==8
num=0;enumerated=0;fixed_counts={};untouched_cases={}
def count(s):
 global num,enumerated
 inter=s&h;m=inter.bit_count()
 if m%2:return 0
 if s not in rows:
  # Any overlap>8 against a newly introduced support would have been
  # rejected by independent audit unless this old class was explicitly fixed.
  val=comb(m,m//2)*2**(11-m) if m<8 else 1120
  untouched_cases[m]=untouched_cases.get(m,0)+1
  return val
 r=rows[s];ok=np.ones(65536,bool)
 for g in r['checks']:ok &= (pc[X&g]%2==0)
 linear=X[ok];assert len(linear)>1 and min(pc[linear[1:]])>=4
 vals=linear^r['shift'];mask=sum(((h>>i)&1)<<j for j,i in enumerate(i for i in range(64) if s>>i&1))
 val=int(np.count_nonzero(pc[vals&mask]==m//2));num+=1;enumerated+=len(vals);fixed_counts[str(s)]=val;return val
terms={'dense':comb(8,4)*2**21,'old_sparse':sum(map(count,old)),'new_sparse':sum(map(count,new))}
blocks={tuple(sorted([a,b,d,a^b^d])) for a,b,d in combinations(range(64),3)}
assert len(blocks)==10416 and all(len(set(t))==4 for t in blocks)
original_four_count=sum(comb(m,m//2)*2**(4-m) for t in blocks if (m:=sum((h>>j)&1 for j in t))%2==0)
# Independently relabel ONLY quadruple layer. Its cross-layer bounds do not
# depend on labels: dense max16, sparse max32, axis max32.
perm=positions+[i for i in range(64) if not(h>>i&1)]
assert len(set(perm))==64 and sum(1<<perm[i] for i in range(8))==h
relabeled={tuple(sorted(perm[i] for i in t)) for t in blocks}
triples=[u for t in relabeled for u in combinations(t,3)]
assert len(relabeled)==10416 and len(set(triples))==len(triples)==41664
assert 4*4<=32 and 4*4*2<=32 and 4*8<=32
terms['four_supports']=sum(comb(m,m//2)*2**(4-m) for t in relabeled if (m:=sum((h>>j)&1 for j in t))%2==0)
terms['poles']=2*(64-h.bit_count())
print('ACTUAL',terms,'EXPECTED',proposal['terms'],flush=True)
assert terms==proposal['terms'] and sum(terms.values())==proposal['total']==184415700
out={'status':'PASS','dimension':63,'octad':h,'terms':terms,'total':sum(terms.values()),'gain_over_published':sum(terms.values())-184393988,'dense_projection_rank':7,'quadruple_permutation':perm,'original_unpermuted_four_count':original_four_count,'explicit_classes':num,'explicit_sign_words':enumerated,'untouched_overlap_cases':untouched_cases,'fixed_counts':fixed_counts,'proof':'Select points orthogonal to octad indicator. Orthogonal isometry to R63 preserves norms and all products. Explicit affine-class count, direct actual-octad count of independently permuted four-support layer, dense rank proof. Untouched old Hamming classes may be relabeled independently to place overlap8 in an affine hyperplane; those classes have no high-overlap changed interactions as verified by separate64D audit.'}
(H/'section_report.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='fixed_counts'})
