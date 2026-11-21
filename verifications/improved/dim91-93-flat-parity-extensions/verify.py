"""Exact replay of the six-/seven-row replacement layer in 91–93D."""
from pathlib import Path
import json, subprocess, sys
P=Path(__file__).resolve().parent; V=P.parent
sparse=V/'dim89-two-color-latitudes'/'abbreviated_parent'
extra=V/'dim88-89-golay-variable-punctures'/'patterns_extra.json'
args=[sparse/'patterns.json',sparse/'weight10.json',extra]
flat=[json.loads(r) for r in subprocess.check_output(['node',str(P/'flat67.cjs'),*map(str,args)],text=True).splitlines()]
dirs=['dim91-z8-asymmetric-planes','dim92-z8-variable-radius-heads-continuation','dim93-z8-variable-radius-heads-continuation']
expected=json.loads((P/'claims.json').read_text())
for d,folder,e,f in zip(range(91,94),dirs,expected,flat,strict=True):
 b=json.loads(subprocess.check_output([sys.executable,str(V/folder/'verify_flat.py')],text=True))
 assert b['dimension']==f['dimension']==d==e['dimension']
 assert b['flat']==323584 and b['old_sparse']==145004544
 assert f['flat']==2613248
 assert f['lowerBound']==e['lower_bound']==b['lower_bound']-b['flat']+f['flat']
 assert e['previous']==b['lower_bound'] and e['gain']==f['gain']==2289664
 print(d,e['lower_bound'],e['gain'])
