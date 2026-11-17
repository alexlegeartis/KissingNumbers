from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P.parent/'dim93-96-e7-rs-root-augmentation/verify.py')],check=True)
for f in ['e7_to_e6_94_20261117/verify.py','e7_mixed_sections95_20261117/verify.py','check.py']:
 subprocess.run([sys.executable,str(P/'e6e7_9495_independent_20261117'/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'94':66573830860,'95':109523504196}
print('PASS: approved 94D and 95D E6/E7 head-exchange bounds')
