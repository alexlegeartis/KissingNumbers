from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P.parent/'dim90-92-hexagonal-phase-layers/verify.py')],check=True)
for f in ['latitude_sparse14_20261117/verify.py','latitude_grid92_20261117/verify.py','latitude14_independent_20261117/check.py','latitude_linear_phases_20261117/verify.py','linear_phase_independent_20261117/check.py','linear_phase_independent_20261117/check_linear.py']:
 subprocess.run([sys.executable,str(P/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'91':21943640934+2*31680*6561+2,'92':21943640934+4*31680*6561+4*6930*2187+6}
print('PASS: direct91D/92D linear-phase latitude layers')
