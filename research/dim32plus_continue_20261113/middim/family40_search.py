"""Whole-family 40D search replacing the entire norm-six core."""
from pathlib import Path
import numpy as np,json,subprocess
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).parent
P=ROOT/'verifications/improved/dim41-46-p48-heads/data'
g=np.load(P/'e8_geometry.npz');mu=g['mu']//4;tiny=np.load(P/'population_table.npz')['tiny'];vertices=[]
for r in [16,14,12,10]:
 for c in range(6561):
  if mu[c]%6!=(18-r)%6 or (r==12 and c==0):continue
  w=116235-40*int(tiny[c]) if r==16 else 9180 if r==14 else 495 if r==12 else int(tiny[c])
  if w:vertices.append((c,r,w))
source=HERE/'family40.txt'
with source.open('w') as f:
 f.write(str(len(vertices))+'\n');f.write(' '.join(map(str,g['powers']))+'\n');f.write(' '.join(map(str,mu))+'\n')
 for row in g['digits']:f.write(' '.join(map(str,row))+'\n')
 for row in vertices:f.write(' '.join(map(str,row))+'\n')
subprocess.run([str(HERE/'family40_search'),str(source),str(HERE/'family40_best.json'),'90'],check=True)
r=json.loads((HERE/'family40_best.json').read_text());chosen=[vertices[i] for i in r['vertices']]
for i,(c,rr,w) in enumerate(chosen):
 for cc,ss,ww in chosen[:i]:
  delta=((g['digits'][c]-g['digits'][cc])%3)@g['powers'];b=rr+ss-18+mu[delta];assert b<=0 or b*b<=rr*ss
assert sum(v[2] for v in chosen)==r['value']
r.update(families=chosen,status='Exact universal finite-family angles and count; inherited populations require replay before a record claim')
(HERE/'family40_result.json').write_text(json.dumps(r,indent=2));print(json.dumps(dict(value=r['value'],families=len(chosen))),flush=True)
