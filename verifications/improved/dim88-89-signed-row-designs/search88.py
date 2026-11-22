from supports import hanani, double
from itertools import combinations
from pathlib import Path
from random import Random
S=[x for x in double(hanani(8),22) if x[-1]<42]
old=set(S)
pool=list(combinations(range(42),5))
rng=Random(210211)
rng.shuffle(pool)
seen=set()
result=[]
for b in pool:
    triples=list(combinations(b,3))
    if any(t in seen for t in triples):continue
    if any(q in old for q in combinations(b,4)):continue
    result.append(b)
    seen.update(triples)
pairs={p for b in result for p in combinations(b,2)}
print(len(result),len(pairs),16*len(result)-4*len(pairs))
Path(__file__).with_name('five.txt').write_text(''.join(' '.join(map(str,b))+'\n' for b in sorted(result)))
