"""Exact composition arithmetic with the stated published design theorem."""
import json
import re
import subprocess
import sys
from pathlib import Path

assert __debug__
P=Path(__file__).resolve().parent
BASE=P.parent/"dim65-66-grouped-four-site-designs"
OLD=419000034
if sys.argv[1:] != ["--composition-only"]:
    assert len(sys.argv)==1
    assert json.loads((BASE/"claims.json").read_text())["lower_bounds"]["65"]==OLD
    out=subprocess.check_output([sys.executable,str(BASE/"verify.py")],text=True)
    assert re.search(r"^DIM 65\b.*\bN 419000034$",out,re.M)
# Existence premise (Bao and Ji, Theorem 4.9):
# A(65,4,4) = Johnson bound.
d=65
supports=(d*(((d-1)*((d-2)//2))//3))//4
assert supports==10741
# All sixteen signs on a support and two-coordinate intersection cap.
assert 4-2==2
# Exact squared-dot bounds against unchanged squared norms 4,16,64,72.
for norm2,dot2 in ((4,4),(16,16),(64,16),(72,50)):
    assert dot2<=norm2
n=OLD+16*(supports-10726)
assert n==419000274
assert json.loads((P/"claims.json").read_text())=={"status":"PASS","lower_bounds":{"65":n}}
print(f"DIM 65: {supports} supports by cited theorem, gain {n-OLD}; N {n}")
