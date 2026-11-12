"""Fresh finite sextet sweeps from November 11 PackingStar seeds."""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from source_paths import REPO, PACKINGSTAR_DEFAULT, packingstar_file
sys.path.insert(0,str(REPO/'verifications/improved/dim27-31-joint-refills'))
from search_sextets import frames,transform,refill,save
p=argparse.ArgumentParser();p.add_argument('--dim',type=int,required=True);p.add_argument('--max-blockers',type=int,default=32);p.add_argument('--frames',type=int,default=8);p.add_argument('--packingstar',type=Path,default=PACKINGSTAR_DEFAULT);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
root=args.output;root.mkdir(parents=True,exist_ok=True)
d=args.dim;t0=time.time();f=np.load(packingstar_file(args.packingstar,d),allow_pickle=False);den=int(f['coordinate_denominator']);x=(f['points']*(2//den)).astype(np.int16)
base=len(x);records=[]
for i,blocks in enumerate(frames()[:args.frames]):
 y=transform(x,blocks)
 if y is None:
  r={'frame':i+1,'skipped':'fractional image'}
 else:
  x,r=refill(x,y,args.max_blockers)
  r.update(frame=i+1,total=len(x),seconds=time.time()-t0)
 records.append(r);print(json.dumps({'dim':d,**r}),flush=True)
 if len(x)>base:
  save(root/f'sextet{d}_{len(x)}.npz',x)
 (root/f'sextet{d}.json').write_text(json.dumps({'dim':d,'baseline':base,'current':len(x),'cutoff':args.max_blockers,'scope':'candidate pending full geometry; negative results restricted to cutoff','frames':records},indent=2)+'\n')
print(json.dumps({'dim':d,'final':len(x),'gain':len(x)-base,'seconds':time.time()-t0}),flush=True)
