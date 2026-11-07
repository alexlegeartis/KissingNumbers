"""Exact small-source upper bound via explicit disjoint conflict-clique covers.
This bounds only the 25 incumbent exceptional supports, not the 171-support pool.
"""
import json,itertools,collections
from pathlib import Path
D=Path(__file__).resolve().parent
c=json.loads((D/'certificate32_352832.json').read_text());ss=sorted({m for m,s in c['exceptional_antipodal_lines']})
triangles=[[0,12,19],[1,16,17],[6,7,8],[2,20,21],[9,18,24],[13,14,15],[3,11,23],[4,10,22]];single=5
assert sorted([single]+sum(triangles,[]))==list(range(25))
def expand(m):
 pos=[j for j in range(32) if m>>j&1]
 return [sum(1<<j for i,j in enumerate(pos) if s>>i&1) for s in range(1<<len(pos))]
def vec(m,s):return [(-2 if s>>j&1 else 2) if m>>j&1 else 0 for j in range(32)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
covers=0
for tri in triangles:
 for i,j in itertools.combinations(tri,2):
  m,k=ss[i],ss[j];common=m&k;a=m^common;b=k^common
  assert common.bit_count() in (5,6,7)
  af=a&-a;bf=b&-b;aa=expand(a^af);bb=expand(b^bf);assert len(aa)==len(bb)
  seen=set();cliques=[]
  for z in expand(common):
   for x,y in zip(aa,bb):
    q=[(m,z|x),(m,z|x|af),(k,z|y),(k,z|y|bf)]
    assert len(set(q))==4 and not seen.intersection(q);seen.update(q)
    vv=[vec(*p) for p in q]
    assert all(sum(t*t for t in v)==32 for v in vv)
    assert all(dot(v,w)>16 for v,w in itertools.combinations(vv,2))
    cliques.append(q)
  assert len(cliques)==128 and seen=={(s,t) for s in [m,k] for t in expand(s)}
  covers+=1
m=ss[single];bit=m&-m;edges=[[(m,s),(m,s|bit)] for s in expand(m^bit)]
assert len(edges)==128 and all(dot(vec(*a),vec(*b))>16 for a,b in edges)
counts=collections.Counter(m for m,s in c['exceptional_antipodal_lines'])
assert sum(2*counts[ss[i]] for t in triangles for i in t)+2*counts[m]==1664
for t in triangles:assert sum(2*counts[ss[i]] for i in t)==192
assert 2*counts[m]==128
out={'all_checks_pass':True,'scope':'Only all 256 sign patterns on each of the same 25 exceptional supports, with other layers fixed','triangles':[[ss[i] for i in t] for t in triangles],'isolated_support':m,'explicit_pair_clique_covers_checked':covers,'each_pair_cover_cliques':128,'triangle_bound':192,'isolated_support_bound':128,'upper_bound':1664,'achieved':1664,'total_configuration':352832}
(D/'pool_bound.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
