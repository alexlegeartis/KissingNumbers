import hashlib,json,subprocess,tempfile
from pathlib import Path
assert __debug__
P=Path(__file__).resolve().parent
W=P/'witness.txt'
assert hashlib.sha256(W.read_bytes()).hexdigest()=='5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f'
with tempfile.TemporaryDirectory() as t:
    exe=str(Path(t)/'verify')
    subprocess.run(['g++','-O3','-std=c++17','-march=native',str(P/'verify.cpp'),'-o',exe],check=True)
    r=subprocess.check_output([exe,str(W)],text=True)
assert len(r.splitlines())==4 and all(l.startswith('PASS') for l in r.splitlines()[:2])
bounds={d:int(x) for d,x in (l.split() for l in r.splitlines()[2:])}
assert bounds==json.loads((P/'claims.json').read_text())['lower_bounds']
print(r,end='')
