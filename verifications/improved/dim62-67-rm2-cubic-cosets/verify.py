"""Replay the explicit cubic sign construction and inherited sparse-layer transfer."""
import os,sys,subprocess,tempfile,json,re
from pathlib import Path
H=Path(__file__).resolve().parent
root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else H.parent/'dim61-67-compatible-sign-families'
assert __debug__
_,a,b=(H/'cubic80.wins').read_text().splitlines()
offsets=list(map(int,a.split()))
group=list(map(int,b.split()))
assert len(offsets)==10 and len(group)==8 and 0 in group
assert {x^y for x in group for y in group}==set(group)
labels={x^y for x in offsets for y in group}
assert len(labels)==80 and labels==set(map(int,(H/'witness.txt').read_text().split()))
with tempfile.TemporaryDirectory() as t:
 for name in ('verify80.cpp','verify_independent.cpp'):
  dst=Path(t)/'verify'
  subprocess.run(['g++','-O3','-std=c++17',str(H/name),'-o',str(dst)],check=True)
  subprocess.run([str(dst),str(H/'witness.txt')],check=True)
result=subprocess.run([sys.executable,str(H/'transfer.py'),str(root)],check=True,capture_output=True,text=True)

actual=json.loads(result.stdout)
section=subprocess.run([sys.executable,str(H/'verify63.py'),str(root.parent),str(H)],check=True,capture_output=True,text=True)
section=json.loads(section.stdout)
assert section==json.loads((H/'verify63.json').read_text())
actual['lower_bounds']['63']=section['total']
actual['lower_bounds']=dict(sorted(actual['lower_bounds'].items()))
assert actual==json.loads((H/'transfer.json').read_text())
print(json.dumps(actual,indent=2))
