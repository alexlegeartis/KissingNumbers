from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
for f in ['three_basis_z8_93_96_20261117/verify.py','three_basis_independent_20261117/check.py','three_basis_independent_20261117/replay_parent.py']:
 subprocess.run([sys.executable,str(P/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'93':40*8**10,'94':60*8**10,'95':84*8**10,'96':3*8**12}
print('PASS: approved93–96 three-crosspolytope-alphabet bounds and parent minimum replay')
