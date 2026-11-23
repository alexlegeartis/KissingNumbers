import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

assert __debug__
P = Path(__file__).resolve().parent
old = P / "old_cubic.txt"
assert hashlib.sha256(old.read_bytes()).hexdigest() == "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f"
with tempfile.TemporaryDirectory() as tmp:
    run = str(Path(tmp) / "verify")
    subprocess.run(["g++", "-O3", "-std=c++17", "-march=native", str(P / "verify.cpp"), "-o", run], check=True)
    log = subprocess.check_output([run], cwd=P, text=True)
rows = log.splitlines()
assert len(rows) == 5 and rows[0].startswith("PASS ")
out = {d: int(n) for d, n in (line.split() for line in rows[1:])}
assert out == json.loads((P / "claims.json").read_text())["lower_bounds"]
print(log, end="")
