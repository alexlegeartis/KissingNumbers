"""Require rejection of four mathematical certificate corruptions."""
from pathlib import Path
import tempfile,subprocess,sys,json,hashlib,numpy as np
P=Path(__file__).resolve().parent
if not __debug__:raise SystemExit('Run without -O')
base=json.loads((P/'certificate.json').read_text());z=np.load(P/'new_points.npz');results=[]
for mode in ['duplicate_point','wrong_head_norm','missing_cap_owner','wrong_total']:
 with tempfile.TemporaryDirectory(prefix='hamming38_negative_') as d:
  d=Path(d);c=json.loads(json.dumps(base));arrays={k:z[k].copy() for k in z.files}
  if mode=='duplicate_point':
   for k in arrays:arrays[k][1]=arrays[k][0]
  elif mode=='wrong_head_norm':arrays['head_integer'][0,0]=2
  elif mode=='missing_cap_owner':c['removed_cap_body_lines'].remove(c['body_line_indices'][0])
  elif mode=='wrong_total':c['claimed_size']+=1
  np.savez_compressed(d/'points.npz',**arrays);c['new_points_sha256']=hashlib.sha256((d/'points.npz').read_bytes()).hexdigest();(d/'certificate.json').write_text(json.dumps(c));r=subprocess.run([sys.executable,'-B',str(P/'verify.py'),'--certificate',str(d/'certificate.json'),'--points',str(d/'points.npz')],capture_output=True,text=True);assert r.returncode!=0,mode;results.append({'case':mode,'rejected':True})
(P/'negative_control_results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
