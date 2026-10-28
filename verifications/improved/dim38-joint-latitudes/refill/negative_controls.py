"""Require both exact proof representations to reject four corrupted witnesses."""
from copy import deepcopy
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent

def main():
    original=json.loads((HERE/'certificate.json').read_text())
    cases={}
    c=deepcopy(original)
    c['doubled_head_lines'][1]=c['doubled_head_lines'][0]
    c['body_line_indices'][1]=c['body_line_indices'][0]
    cases['duplicate grid family']=c
    c=deepcopy(original)
    c['additional_latitude_points'][0]['body_sign'] *= -1
    cases['incorrect cube body sign']=c
    c=deepcopy(original)
    c['removed_cap_body_lines'].pop()
    c['claimed_size'] += 6
    cases['undeleted same-body caps']=c
    c=deepcopy(original)
    c['claimed_size'] += 1
    cases['false claimed count']=c
    environment=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
    with tempfile.TemporaryDirectory(prefix='kissing38_negative_') as temporary:
        path=Path(temporary)/'bad.json'
        for name,c in cases.items():
            path.write_text(json.dumps(c))
            for script in ('verify.py','verify_coordinates.py'):
                run=subprocess.run([sys.executable,'-B',str(HERE/script),'--certificate',str(path)],
                    text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,env=environment)
                assert run.returncode != 0 and 'ValueError:' in run.stdout, (name,script,run.stdout)
                print('REJECTED:',name,'by',script,flush=True)
    print('PASS: all eight independent rejection checks.')

if __name__=='__main__':
    main()
