#!/usr/bin/env python3
"""Recompute all counts and exact modified pair-category certificates."""
import argparse,json,subprocess,sys
from pathlib import Path
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--siblings',type=Path,default=P.parent)
args=parser.parse_args()
if not __debug__:raise SystemExit('Run with assertions enabled.')
claims=json.loads((P/'claims.json').read_text())
got={}
dims={'d54':54,'d55':55,'d57_59':57,'d60':60,'d57trine':57,'d55trine':55,'d59trine':59,'d60cube':60}
for p in dims:
    result=subprocess.run([sys.executable,str(P/p/'verify.py'),'--siblings',
                           str(args.siblings.resolve())],check=True,capture_output=True,text=True)
    print(result.stdout,end='')
    for line in result.stdout.splitlines():
        q=json.loads(line)
        d=str(q.get('dimension',dims[p]))
        got[d]=dict(previous=q['previous'], total=q['lower_bound'])
assert got==claims,(got,claims)
print('Exact geometry replay and all five registry counts confirmed.')
