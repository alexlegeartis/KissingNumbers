#!/usr/bin/env python3
"""Check both signed predecessors and add the new norm-four tails."""
from pathlib import Path
import json,runpy,sys
from hashlib import sha256
from norm4 import norm4
assert __debug__
P=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
verify85=BASE/'dim85-gamma72-signed-apron'
assert all(sha256((verify85/k).read_bytes()).hexdigest()==v
           for k,v in json.loads((P/'pins85.json').read_text()).items())
answers=json.loads((P/'claims.json').read_text())['lower_bounds']
for n,file,origin,bf in [
 (84,P/'signed84.py',P,'apron84.json'),
 (85,verify85/'verify.py',verify85,'claims.json')]:
    sys.argv=[str(file),str(BASE)]
    sys.path.insert(0,str(origin))
    K=runpy.run_path(sys.argv[0])
    sys.path.remove(str(origin))
    p=json.loads((origin/bf).read_text())['lower_bounds'][str(n)]
    extra=norm4(n,K,origin)
    assert answers[str(n)]==p+extra
    print('PASS signed norm-four tails',n,p+extra,'increment',extra,flush=True)
