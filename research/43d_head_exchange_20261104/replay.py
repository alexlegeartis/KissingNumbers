"""Replay an exact 43D three-for-five exchange without modifying the repository."""
import base64,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
if not __debug__:raise SystemExit('Run without -O; assertions are part of the verifier.')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PACKAGE=ROOT/'verifications/improved/dim41-46-p48-heads';SRC=PACKAGE/'partial43';CHECK=PACKAGE/'continuation43/verify_coordinates.py';data=json.loads((HERE/'delta.json').read_text());baseline_bytes=base64.b64decode((HERE/'baseline_2919782.npz.b64').read_bytes())
assert hashlib.sha256(baseline_bytes).hexdigest()=='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4'
with tempfile.TemporaryDirectory(prefix='k43-exchange-') as temp:
 root=Path(temp);dest=root/'partial43';shutil.copytree(SRC,dest);base=root/'baseline_2919782.npz';base.write_bytes(baseline_bytes);old=dict(np.load(base,allow_pickle=False));deleted=np.array(data['removed_indices'],np.int64);heads=np.array(data['head_numerators'],np.int64);bodies=np.array(data['added_body_indices'],np.int64);D=int(old['head_denominator']);assert D==data['head_denominator'] and len(bodies)==5 and len(deleted)==3 and np.array_equal(deleted,np.unique(deleted))
 new=dict(old);new['body_indices']=np.r_[np.delete(old['body_indices'],deleted),bodies];new['head_numerators']=np.vstack((np.delete(old['head_numerators'],deleted,axis=0),heads));points=dest/'explicit_points_43.npz';np.savez_compressed(points,**new);delta=root/'delta.npz';np.savez_compressed(delta,deleted_indices=deleted,added_body_indices=bodies,added_head_numerators=heads,head_denominator=np.int64(D))
 env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1';checks=[[sys.executable,str(dest/'verify.py')],[sys.executable,str(CHECK),'--root',str(root),'--points',str(points),'--previous',str(base),'--delta',str(delta)]]
 for cmd in checks:
  result=subprocess.run(cmd,check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env);report=json.loads(result.stdout);assert report['status']=='PASS' and report['points']==2919784 and report['explicit_tiny_points']==3449;print(json.dumps(report,indent=2))
 print('PASS: both exact replays certify the 2,919,784-point construction, subject to the inherited parent and population inputs.')
