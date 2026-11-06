from pathlib import Path
from decimal import Decimal
import numpy as np,json,hashlib
if not __debug__:raise SystemExit('Run without -O: exact verification uses assertions.')
D=Path(__file__).parent;p=D/'source/coordinates841.txt';rows=[[Decimal(s) for s in line.split()] for line in p.read_text().splitlines() if line.strip()];assert len(rows)==841 and all(len(row)==12 for row in rows)
scale=10**max(-x.as_tuple().exponent for row in rows for x in row)
X=np.array([[int(x.as_integer_ratio()[0]*(scale//x.as_integer_ratio()[1])) for x in row] for row in rows],dtype=object);n=(X*X).sum(1);assert all(v>0 for v in n)
G=X@X.T;bad=(G>0)&(4*G*G>n[:,None]*n);np.fill_diagonal(bad,False);assert not bad.any()
r={'points':841,'dimension':12,'pairs_checked':841*840//2,'violations':0,'arithmetic':'Python arbitrary-precision integers on exact printed decimal rays','sha256':hashlib.sha256(p.read_bytes()).hexdigest()};(D/'pole_verification.json').write_text(json.dumps(r,indent=2));print(r)
