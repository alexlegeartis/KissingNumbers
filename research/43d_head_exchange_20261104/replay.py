"""Replay a finite43D exchange without modifying the registered certificate."""
import argparse,base64,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
if not __debug__:raise SystemExit('Run without -O; assertions are part of the verifier.')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--delta',type=Path,default=HERE/'delta.json');args=parser.parse_args()
PACKAGE=ROOT/'verifications/improved/dim41-46-p48-heads';SRC=PACKAGE/'partial43';CHECK=PACKAGE/'continuation43/verify_coordinates.py'
data=json.loads(args.delta.read_text());baseline_bytes=base64.b64decode((HERE/'baseline_2919782.npz.b64').read_bytes())
assert hashlib.sha256(baseline_bytes).hexdigest()=='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4'
with tempfile.TemporaryDirectory(prefix='k43-exchange-') as temp:
 root=Path(temp);dest=root/'partial43';shutil.copytree(SRC,dest)
 base=root/'baseline_2919782.npz';base.write_bytes(baseline_bytes);old=dict(np.load(base,allow_pickle=False))
 if 'removed_indices' in data:
  deleted=np.asarray(data['removed_indices'],np.int64);bodies=np.asarray(data['added_body_indices'],np.int64)
 else:
  deleted=np.asarray([data['remove']],np.int64);bodies=np.asarray([x[0] for x in data['add']],np.int64)
 heads=np.asarray(data.get('head_numerators',data.get('added_head_numerators')),np.int64);D=int(old['head_denominator'])
 assert D==data['head_denominator'] and len(heads)==len(bodies) and len(bodies)>len(deleted)
 assert np.array_equal(deleted,np.unique(deleted)) and np.all((deleted>=0)&(deleted<len(old['body_indices'])))
 new=dict(old);new['body_indices']=np.r_[np.delete(old['body_indices'],deleted),bodies]
 new['head_numerators']=np.vstack((np.delete(old['head_numerators'],deleted,axis=0),heads))
 points=dest/'explicit_points_43.npz';np.savez_compressed(points,**new)
 delta=root/'delta.npz';np.savez_compressed(delta,deleted_indices=deleted,added_body_indices=bodies,added_head_numerators=heads,head_denominator=np.int64(D))
 model=np.load(dest/'configuration_43.npz');total=sum(map(int,model['weights']))+len(model['pure_directions'])+len(new['body_indices'])
 env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1'
 checks=[[sys.executable,str(dest/'verify.py')],[sys.executable,str(CHECK),'--root',str(root),'--points',str(points),'--previous',str(base),'--delta',str(delta)]]
 for cmd in checks:
  result=subprocess.run(cmd,check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
  report=json.loads(result.stdout);assert report['status']=='PASS' and report['points']==total and report['explicit_tiny_points']==len(new['body_indices'])
  print(json.dumps(report,indent=2))
 print(f'PASS: both exact replays certify the {total:,}-point construction, subject to inherited parent and population inputs.')
