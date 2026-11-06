"""Replay two product-level rejection controls against the independent C++ code."""
from pathlib import Path
import json,subprocess,tempfile
H=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
c=json.loads((H.parent/'orbit_certificate.json').read_text());v=c['representatives'][0];L=c['norm_lower'];results={}
with tempfile.TemporaryDirectory(prefix='kissing39_controls_') as directory:
 exe=Path(directory)/'relative'
 subprocess.run(['g++','-O3','-std=c++17',str(H/'relative.cpp'),'-o',str(exe)],check=True)
 for name,rows,group in [('duplicate_representative',[v,v],[0]),('invalid_self_orbit',[v],[0,1])]:
  data=f'{len(rows)} 39 {len(group)} {L}\n'+' '.join(map(str,group))+'\n'+'\n'.join(' '.join(map(str,row)) for row in rows)+'\n'
  p=subprocess.run([str(exe)],input=data,text=True,capture_output=True)
  assert p.returncode==1 and 'VIOLATION' in p.stderr,(name,p.returncode,p.stderr)
  results[name]={'rejected':True,'stderr':p.stderr.strip()}
(H/'negative_controls.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS: both independent product controls rejected')
