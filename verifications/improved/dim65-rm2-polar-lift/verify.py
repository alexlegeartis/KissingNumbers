#!/usr/bin/env python3
import hashlib,json,subprocess,tempfile
from pathlib import Path
assert __debug__
H = Path(__file__).resolve().parent
P = H / 'witness.txt'
assert hashlib.sha256(P.read_bytes()).hexdigest() == '5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f'
with tempfile.TemporaryDirectory() as t:
    f = str(Path(t) / 'lift65')
    subprocess.run(['g++','-O3','-std=c++17','-march=native',str(H/'lift65.cpp'),'-o',f],check=True)
    line = subprocess.check_output([f,str(P)],text=True).strip()
    assert int(line.split()[-1]) == json.loads((H/'claims.json').read_text())['lower_bound']
    print(line)
