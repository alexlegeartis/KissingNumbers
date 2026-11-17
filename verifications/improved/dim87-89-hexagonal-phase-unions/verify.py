from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
Q=P.parent/'dim68-91-lattice-code-continuation'
for f in ['ternary_qr48_independent_20261117/replay_generator.py','ternary_head88_independent_20261117/verify.py','ternary_heads87_89_independent_20261117/verify.py']:
 subprocess.run([sys.executable,str(Q/f)],check=True,cwd=Q)
subprocess.run([sys.executable,str(P.parent/'dim90-92-hexagonal-phase-layers/verify.py')],check=True)
for f in ['hexagon_87_89_consequence_20261117/verify.py','hexagon_87_89_independent_20261117/audit.py']:
 subprocess.run([sys.executable,str(P/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'87':7266268830,'88':9590791776,'89':12690155704}
print('PASS: 87–89D dense, nonlinear-phase and pure-head unions')
