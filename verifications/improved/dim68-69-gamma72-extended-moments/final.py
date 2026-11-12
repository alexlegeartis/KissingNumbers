#!/usr/bin/env python3
"""Register both bounds only after the local exact certificate check passes."""
from pathlib import Path
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def results():
    replay=subprocess.run([sys.executable,str(HERE/'verify.py')],
                          check=True,capture_output=True,text=True)
    reports=json.loads(replay.stdout)
    return [dict(dimension=r['dimension'],total=r['bound'],previous=r['prior_bound'],gain=r['gain'])
            for r in reports]

if __name__=='__main__':
    print(json.dumps(results(),indent=2))
