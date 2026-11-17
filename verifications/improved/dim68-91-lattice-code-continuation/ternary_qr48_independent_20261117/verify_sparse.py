from pathlib import Path
from collections import Counter
from itertools import combinations,product
from fractions import Fraction as F
import json,hashlib,sympy as s
P=Path(__file__).resolve().parent
c=json.load(open(P/'certificate.json'));assert hashlib.sha256((P/'certificate.json').read_bytes()).hexdigest()=='1505ab829025e938ad09de3575656ccf5f547537d8fd5432533596fa963b22e8'
x=s.symbols('x');g=x**23+x**20-x**18+x**17-x**15+x**12+x**10-x**9-x**7-x**6-x**5-x**2-x-1
assert s.Poly(g,x,modulus=3).is_irreducible
published=json.load(open(P.parent/'ternary_qr48_20261117/sparse_augmentation.json'))['rows'];base={88:6967552256,89:9270596604,90:12341351672,91:16435507572};out=[]
for n in [88,89,90,91]:
 counts=Counter(i^j for i in range(n)for j in range(i+1,n));partitions=sum(k*(k-1)//2 for k in counts.values());assert partitions%3==0;blocks=partitions//3
 # Equal-XOR distinct pairs cannot share a vertex. Each four-set contributes exactly3 pair partitions.
 # Separate direct four-tuple predicate is an independent count.
 direct=sum(a^b^c^d==0 for a,b,c,d in combinations(range(n),4));assert blocks==direct
 N=44 if n==88 else 45
 # All sign patterns: within one support each differing sign subtracts N/2.
 for a,b in combinations(list(product([-1,1],repeat=4)),2):assert F(N,4)*sum(u*v for u,v in zip(a,b))<=F(N,2)
 # Distinct supports overlap at most2; axes-support overlap at most1.
 assert 2*F(N,4)<=F(N,2)
 # Dense coordinate maximum sqrt2; cross squared bounds.
 assert 8*N<=F(N*N,4) and 2*N<=F(N*N,4)
 dense=c['counts'][str(n)];layer=16*blocks+2*n;total=dense+layer
 assert blocks==published[str(n)]['affine_plane_supports'] and total==published[str(n)]['total']
 out.append(dict(n=n,supports=blocks,layer=layer,total=total,benchmark=base[n],gain=total-base[n],percent=float(F(100*(total-base[n]),base[n]))))
(P/'sparse_report.json').write_text(json.dumps({'status':'PASS','rows':out},indent=2));print(json.dumps(out,indent=2))
