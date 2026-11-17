import itertools,json,hashlib,functools,operator
from pathlib import Path
p=Path(__file__).parent;d=json.loads((p/'input.json').read_text())
xor=lambda x:functools.reduce(operator.xor,x,0)
Q=list(map(tuple,d['group_quads'])); assert len(Q)==len(set(Q))==30; assert all(len(set(q))==4 and all(0<=i<10 for i in q)for q in Q); assert sorted(t for q in Q for t in itertools.combinations(sorted(q),3))==list(itertools.combinations(range(10),3))
rounds=[list(map(tuple,r))for r in d['near_perfect_matchings']];edges=[tuple(sorted(e))for r in rounds for e in r];assert len(set(edges))==36 and set(edges)==set(itertools.combinations(range(9),2))
for r in rounds:assert len(set(sum((list(e)for e in r),[])))==8
words=[w for w in itertools.product(range(4),repeat=4)if xor(w)==0];assert set(words)==set(map(tuple,d['pair_label_code']))
assert all(sum(a!=b for a,b in zip(x,y))>=2 for x,y in itertools.combinations(words,2))
assert all(len(set(x)&set(y))<=2 for x,y in itertools.combinations(Q,2))
signs=d['sign_masks'];assert len(signs)==len(set(signs))==128 and all(0<=x<256 for x in signs);assert min((x^y).bit_count()for x,y in itertools.combinations(signs,2))==2
for row in d['rows']:
 n=row['dimension']; k=min(n-81,9);supports=[]
 for q in Q:
  for r in rounds:
   for w in words:
    s=tuple(sorted(9*g+j for g,label in zip(q,w)for j in r[label]))
    if all(j<81+k for j in s):supports.append(s)
 assert len(supports)==len(set(supports))==row['support8_count']
 assert hashlib.sha256(json.dumps(sorted(supports)).encode()).hexdigest()==row['supports_sha256']
 assert row['total']==row['previous_total']+128*len(supports)
print('PASS independent support reconstruction, matching/label/sign intersection lemmas and all five totals')
