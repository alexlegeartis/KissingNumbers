import json,itertools
from pathlib import Path
p=Path(__file__).parent;d=json.loads((p/'input.json').read_text());S=d['sign_masks'];words=[tuple((x>>i)&1 for i in range(10))for x in S]
assert len(set(words))==72 and all(0<=x<1024 for x in S)
ds=[sum(a!=b for a,b in zip(x,y))for x,y in itertools.combinations(words,2)];assert len(ds)==2556 and min(ds)==3
for row in d['rows']:
 n=row['dimension'];support=min(9,n-81)*59049;dense={87:6973568802,88:9298091736,89:12397455648,90:20920706406,91:20920706406}[n]
 assert row['total']==dense+72*support+2*n*(n-1)
print('PASS: all 2556 distances and five totals; minimum',min(ds),'maximum',max(ds))
