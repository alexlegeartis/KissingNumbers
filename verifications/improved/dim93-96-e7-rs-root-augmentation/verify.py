from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P.parent/'dim93-96-three-crosspolytope-alphabets/verify.py')],check=True)
for f in ['three_basis_final93_96_20261117/e7/verify.py','three_basis_final93_96_20261117/verify_augmentation.py','three_basis_final_independent_20261117/check.py']:
 subprocess.run([sys.executable,str(P/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'93':42951379608,'94':64426347212,'95':107376020548,'96':206160530240}
print('PASS: approved93–96 E7 and RS/root augmentation bounds')
