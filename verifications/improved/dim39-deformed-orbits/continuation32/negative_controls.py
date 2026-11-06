"""Require the primary geometry checker to reject four corrupted certificates."""
import json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
original=json.loads((H/'orbit_certificate.json').read_text());results=[]
base=json.loads((H/'base_certificate.json').read_text());supports=[[i for i,v in enumerate(s) if v=='1'] for s in (H/'source/supports_3424.txt').read_text().splitlines()]
for edit in base['support_replacements']:supports[edit['remove_index']]=edit['add']
for name,expected in [('duplicate_orbit','distinct orbit pair'),('zero_norm','Exact squared norms'),('duplicate_group','Distinct group elements'),('middle_collision','Middle support angle')]:
 with tempfile.TemporaryDirectory(prefix='kissing39_actual_control_') as temp:
  target=Path(temp)/'source';shutil.copytree(H,target,ignore=shutil.ignore_patterns('__pycache__','independent'));c=json.loads(json.dumps(original))
  if name=='duplicate_orbit':c['representatives'][1]=c['representatives'][0].copy()
  elif name=='zero_norm':c['representatives'][0]=[0]*39
  elif name=='duplicate_group':c['group_words'][1]=c['group_words'][0]
  else:
   c['representatives'][0]=[35355339 if j in supports[0] else 0 for j in range(39)]
  (target/'orbit_certificate.json').write_text(json.dumps(c));p=subprocess.run([sys.executable,str(target/'verify.py')],env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True)
  assert p.returncode!=0 and expected in p.stderr,(name,p.stdout,p.stderr)
  results.append({'control':name,'rejected':True,'checked_error':expected});print(name,'rejected',flush=True)
(H/'negative_control_results.json').write_text(json.dumps(results,indent=2)+'\n');print('PASS: four corruptions rejected')
