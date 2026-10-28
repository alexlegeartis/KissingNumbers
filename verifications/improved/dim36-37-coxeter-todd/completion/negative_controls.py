"""Reject a duplicate orbit and a false count in both complete exact proofs."""
from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent


def main():
    assert __debug__, 'Run without -O'
    original=json.loads((HERE/'certificate.json').read_text())
    cases={}
    bad=deepcopy(original);bad['phase_orbits'][1]=bad['phase_orbits'][0]
    cases['duplicate phase orbit']=bad
    bad=deepcopy(original);bad['points_36']+=1;bad['points_37']+=1
    cases['false point counts']=bad
    with tempfile.TemporaryDirectory(prefix='phase_negative_') as directory:
        root=Path(directory);target=root/'completion';target.mkdir()
        (root/'phase').mkdir()
        shutil.copyfile(HERE.parent/'phase/exact_lattice.py',root/'phase/exact_lattice.py')
        for script in ('verify.py','verify_coordinates.py'):
            shutil.copyfile(HERE/script,target/script)
        for label,bad in cases.items():
            (target/'certificate.json').write_text(json.dumps(bad))
            for script in ('verify.py','verify_coordinates.py'):
                result=subprocess.run([sys.executable,'-B',str(target/script)],
                    stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,
                    env=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2'))
                assert result.returncode!=0 and 'AssertionError' in result.stdout,(label,script,result.stdout)
                print('REJECTED:',label,'by',script,flush=True)
    print('PASS: all four exact-proof rejection checks')


if __name__=='__main__':
    main()
