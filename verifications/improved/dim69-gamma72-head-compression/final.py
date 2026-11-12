#!/usr/bin/env python3
"""Register the 69D value only after its exact orthonormal proof passes."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def results():
    spec=importlib.util.spec_from_file_location('gamma72_head69_orthonormal',HERE/'verify_orthonormal.py')
    verifier=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    report=verifier.check(HERE/'orthonormal_certificate.json')
    return [{'dimension':report['dimension'],'total':report['bound'],
             'previous':report['prior_bound'],'gain':report['gain']}]

if __name__=='__main__':
    subprocess.run([sys.executable,str(HERE/'verify.py')],check=True,capture_output=True,text=True)
    print(json.dumps(results(),indent=2))
