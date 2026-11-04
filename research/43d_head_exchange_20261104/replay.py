"""Replay an exact 43D one-for-two exchange without modifying the repository."""
import hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
if not __debug__:raise SystemExit('Run without -O; assertions are part of the verifier.')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PACKAGE=ROOT/'verifications/improved/dim41-46-p48-heads';SRC=PACKAGE/'partial43';CHECK=PACKAGE/'continuation43/verify_coordinates.py';data=json.loads((HERE/'delta.json').read_text());base=SRC/'explicit_points_43.npz'
assert hashlib.sha256(base.read_bytes()).hexdigest()=='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4'
with tempfile.TemporaryDirectory(prefix='k43-exchange-') as temp:
 root=Path(temp);dest=root/'partial43';shutil.copytree(SRC,dest);old=dict(np.load(base,allow_pickle=False));i=data['remove'];heads=np.array(data['head_numerators'],np.int64);bodies=np.array([x[0] for x in data['add']],np.int64);D=int(old['head_denominator']);assert D==data['head_denominator'] and len(bodies)==2
 new=dict(old);new['body_indices']=np.r_[np.delete(old['body_indices'],i),bodies];new['head_numerators']=np.vstack((np.delete(old['head_numerators'],i,axis=0),heads));points=dest/'explicit_points_43.npz';np.savez_compressed(points,**new);delta=root/'delta.npz';np.savez_compressed(delta,deleted_indices=np.array([i],np.int64),added_body_indices=bodies,added_head_numerators=heads,head_denominator=np.int64(D))
 env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1';checks=[[sys.executable,str(dest/'verify.py')],[sys.executable,str(CHECK),'--root',str(root),'--points',str(points),'--previous',str(base),'--delta',str(delta)]]
 for cmd in checks:
  result=subprocess.run(cmd,check=True,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env);report=json.loads(result.stdout);assert report['status']=='PASS' and report['points']==2919783 and report['explicit_tiny_points']==3448;print(json.dumps(report,indent=2))
 print('PASS: both exact replays certify the 2,919,783-point construction, subject to the inherited parent and population inputs.')
