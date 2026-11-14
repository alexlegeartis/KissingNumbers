import argparse,json,itertools,time
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parents[2];rows=[]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
files={26:'dim26-27-radical-orbits/data/coordinates26_201188.npz',27:'dim26-27-radical-orbits/data/coordinates27_209884.npz',28:'dim28-29-november-sextets/data/coordinates28_222120.npz',29:'dim29-30-native-partial-design/data/coordinates29_231814.npz',30:'dim29-30-native-partial-design/data/coordinates30_248180.npz',31:'dim31-annealed-design/data/coordinates31_251420.npz'}
for d,f in files.items():
 p=R/'verifications/improved'/f
 if not p.exists():print('MISSING',p,flush=True);continue
 z=np.load(p);print('INPUT',d,z.files,flush=True)
 A=z['rational'] if 'rational' in z else z['points'];B=z['radical'] if 'radical' in z else np.zeros_like(A)
 best=(0,None)
 for i in range(d):
  n=int(np.sum((A[:,i]==0)&(B[:,i]==0)))
  if n>best[0]:best=(n,[i])
 for i,j in itertools.combinations(range(d),2):
  for s in (-1,1):
   n=int(np.sum((A[:,i]==s*A[:,j])&(B[:,i]==s*B[:,j])))
   if n>best[0]:best=(n,[i,j,s])
 row=dict(dimension=d,target=d-1,count=best[0],normal=best[1]);rows.append(row);print('BEST_SECTION',json.dumps(row),flush=True)
args.out.write_text(json.dumps(rows,indent=2)+'\n')
