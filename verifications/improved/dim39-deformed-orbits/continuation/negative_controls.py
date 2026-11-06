"""Require both arithmetic routes to reject four independently corrupted inputs."""
import json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
H=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
original=json.loads((H/'orbit_certificate.json').read_text());results=[]
for name in ('duplicate_orbit','false_count','zero_norm','duplicate_group'):
 for checker in ('verify.py','verify_direct.py'):
  with tempfile.TemporaryDirectory(prefix='kissing39_negative_') as temp:
   root=Path(temp)/'source';shutil.copytree(H,root,ignore=shutil.ignore_patterns('__pycache__'))
   c=json.loads(json.dumps(original))
   if name=='duplicate_orbit':c['representatives'][1]=c['representatives'][0].copy()
   elif name=='false_count':c['claimed_size']+=1
   elif name=='zero_norm':c['representatives'][0]=[0]*39
   else:c['group_words'][1]=c['group_words'][0]
   (root/'orbit_certificate.json').write_text(json.dumps(c))
   env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
   p=subprocess.run([sys.executable,str(root/checker)],cwd=root,env=env,capture_output=True,text=True)
   assert p.returncode!=0 and 'Traceback' in p.stderr,(name,checker,p.stdout,p.stderr)
   results.append(dict(control=name,checker=checker,rejected=True));print(name,checker,'rejected',flush=True)
(H/'negative_control_results.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS: all eight corruptions rejected')
