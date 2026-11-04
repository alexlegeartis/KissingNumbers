"""Install the independently replayed engine delta, preserving prior artifacts.

Run only after replay.py --delta engine43_combined_delta.json passes both paths.
This updates the current ledger, exact delta and package integrity manifests.
"""
from pathlib import Path
import base64,hashlib,io,json,os,shutil,subprocess,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
PACKAGE=ROOT/'verifications/improved/dim41-46-p48-heads';CONT=PACKAGE/'continuation43'

def write_text(path,text):
 ending=b'\r\n' if path.exists() and b'\r\n' in path.read_bytes() else b'\n'
 path.write_bytes(text.replace('\r\n','\n').replace('\n',ending.decode()).encode())

def write_json(path,value):write_text(path,json.dumps(value,indent=2)+'\n')

def reports(text):
 decoder=json.JSONDecoder();result=[]
 while text.lstrip().startswith('{'):
  obj,n=decoder.raw_decode(text.lstrip());result.append(obj);text=text.lstrip()[n:]
 return result

def main():
 source=ROOT/'research/43d_head_exchange_20261104';raw=base64.b64decode((source/'baseline_2919782.npz.b64').read_bytes())
 if hashlib.sha256(raw).hexdigest()!='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4':raise ValueError('archived predecessor')
 current=PACKAGE/'partial43/explicit_points_43.npz'
 if current.read_bytes()!=raw:raise ValueError('Canonical point file changed; inspect concurrent work before installing')
 # A stored PASS receipt is not authority. Re-run both independent geometric
 # representations before changing any canonical data or registry entry.
 env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1'
 replay=subprocess.run([sys.executable,str(source/'replay.py'),'--delta',str(OUT/'engine43_combined_delta.json')],
  check=True,text=True,capture_output=True,env=env)
 proof=reports(replay.stdout)
 if len(proof)!=2 or any(r.get('status')!='PASS' or r.get('points')!=2919785 for r in proof):raise ValueError('Both completed independent replays are required')
 old=dict(np.load(io.BytesIO(raw),allow_pickle=False));delta=json.loads((OUT/'engine43_combined_delta.json').read_text());D=int(old['head_denominator'])
 deleted=np.array(delta['removed_indices'],np.int64);bodies=np.array(delta['added_body_indices'],np.int64);heads=np.array(delta['head_numerators'],np.int64)
 assert len(deleted)==5 and len(bodies)==8 and D==delta['head_denominator']
 final=dict(old);final['body_indices']=np.r_[np.delete(old['body_indices'],deleted),bodies];final['head_numerators']=np.vstack((np.delete(old['head_numerators'],deleted,axis=0),heads))
 # The generic engine selected the same final point set as the independently
 # developed three-step construction, in a different deterministic order.
 parent=json.loads((source/'delta.json').read_text())
 parent_ix=np.r_[np.delete(old['body_indices'],parent['removed_indices']),parent['added_body_indices']]
 parent_h=np.vstack((np.delete(old['head_numerators'],parent['removed_indices'],axis=0),parent['head_numerators']))
 key=lambda ix,h:{(int(i),*map(int,z)) for i,z in zip(ix,h)}
 assert key(parent_ix,parent_h)==key(final['body_indices'],final['head_numerators'])
 history=CONT/'history_2919782';history.mkdir(exist_ok=True)
 archived={'explicit_points_43.npz':current,'delta.npz':CONT/'delta.npz','README.md':CONT/'README.md',
  'ambient_verification.json':CONT/'ambient_verification.json','coordinate_verification.json':CONT/'coordinate_verification.json',
  'reproduction.json':CONT/'reproduction.json','43_CERTIFICATION.md':PACKAGE/'43_CERTIFICATION.md'}
 for name,path in archived.items():
  dst=history/name
  if dst.exists() and dst.read_bytes()!=path.read_bytes():raise ValueError('Historical archive mismatch')
  if not dst.exists():shutil.copyfile(path,dst)
 # Compose against the portable3408-point predecessor used by the package.
 previous=dict(np.load(CONT/'previous_points.npz',allow_pickle=False));pd=int(previous['head_denominator']);assert D%pd==0
 prior={(int(b),*map(int,z*(D//pd))):i for i,(b,z) in enumerate(zip(previous['body_indices'],previous['head_numerators']))}
 keep=set();extra=[]
 for i,(b,z) in enumerate(zip(final['body_indices'],final['head_numerators'])):
  k=(int(b),*map(int,z))
  if k in prior:keep.add(prior[k])
  else:extra.append(i)
 combined={'deleted_indices':np.array(sorted(set(range(len(previous['body_indices'])))-keep),np.int64),
  'added_body_indices':final['body_indices'][extra],'added_head_numerators':final['head_numerators'][extra],
  'head_denominator':np.int64(D)}
 np.savez_compressed(CONT/'delta.npz',**combined);np.savez_compressed(current,**final)
 assert hashlib.sha256(current.read_bytes()).hexdigest()==proof[1]['points_sha256']
 write_json(PACKAGE/'partial43/replay.json',proof[0]);write_json(CONT/'ambient_verification.json',proof[0]);write_json(CONT/'coordinate_verification.json',proof[1])
 model=np.load(PACKAGE/'partial43/configuration_43.npz');total=int(sum(map(int,model['weights'])))+len(model['pure_directions'])+len(final['body_indices'])
 ledger=json.loads((PACKAGE/'bounds.json').read_text());row=next(r for r in ledger['constructions'] if r['dimension']==43)
 assert row['total']==2919782 and total==2919785
 row.update(total=total,explicit_body_points=len(final['body_indices']),body_points=total-row['pure_points'])
 ledger['checkpoint_date']='2026-11-04';write_json(PACKAGE/'bounds.json',ledger)
 provenance=json.loads((PACKAGE/'UPDATE_PROVENANCE.json').read_text())
 provenance['body_head43_20261104']={'previous_total':2919782,'total':total,'gain':3,'tiny_points':3450,
  'source_branch_tip':'f050d67cb0693d7b0d706b80b71ac91314e5288b','original_discovery':'Separate finite43D searches; see research/43d_head_exchange_20261104',
  'engine':'search.body_head_exchange independently generated all three steps from body/head/envelope inputs, without a supplied delta or target count',
  'engine_delta_sha256':hashlib.sha256((OUT/'engine43_combined_delta.json').read_bytes()).hexdigest(),
  'point_certificate_sha256':hashlib.sha256(current.read_bytes()).hexdigest(),
  'combined_deleted_from3408':len(combined['deleted_indices']),'combined_added_to3408':len(extra),
  'population_scope':'Unchanged inherited complete-family identities and P48n minimum/design inputs; no fresh large-shell enumeration.'}
 write_json(PACKAGE/'UPDATE_PROVENANCE.json',provenance)
 write_json(OUT/'registry43_installation.json',provenance['body_head43_20261104'])
 print(json.dumps(provenance['body_head43_20261104'],indent=2))
if __name__=='__main__':main()
