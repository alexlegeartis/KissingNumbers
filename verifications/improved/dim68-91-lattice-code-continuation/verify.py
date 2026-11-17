"""Run unchanged finite proof checkers; published source premises remain explicit."""
from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
for path in ['ternary_qr48_independent_20261117/replay_generator.py','ternary_qr48_independent_20261117/verify_sparse.py','ternary_head88_independent_20261117/verify.py','ternary_heads87_89_independent_20261117/verify.py','mixed_shell71_20261117/verify71.py']:
    subprocess.run([sys.executable,str(P/path)],cwd=P,check=True)
subprocess.run([sys.executable,str(P/'verify_lattice.py')],cwd=P,check=True)
for path in ['rs9_sparse_independent_20261117/verify.py','rs9_nonlinear72_independent_20261117/verify.py','rs9_sqs10_independent_20261117/verify.py']:
    subprocess.run([sys.executable,str(P/path)],cwd=P,check=True)
expected={'68':472598888,'69':780388888,'71':2623578936,'87':7000788678,'88':9329710944,'89':12433498768,'90':20961198018,'91':20961198378}
assert json.loads((P/'claims.json').read_text())==expected
for row in json.loads((P/'rs9_sqs10_independent_20261117/input.json').read_text())['rows']:
    assert expected[str(row['dimension'])]==row['total']
print('PASS: all eight approved published lower bounds and unchanged independent certificates.')
