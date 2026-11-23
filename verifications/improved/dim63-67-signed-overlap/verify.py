import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

assert __debug__
P = Path(__file__).resolve().parent
assert hashlib.sha256((P / "cubic80.txt").read_bytes()).hexdigest() == "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f"
with tempfile.TemporaryDirectory() as d:
    D = Path(d)
    for name in ("cubic80.txt", "polar.txt"):
        shutil.copyfile(P / name, D / name)
    data = gzip.decompress((P / "signed_witness.txt.gz").read_bytes())
    assert hashlib.sha256(data).hexdigest() == "5173fabd124bd23e077d00c956534e7f56543101b28903cbac017ed252c34be7"
    (D / "signed_witness.txt").write_bytes(data)
    binary = str(D / "check")
    subprocess.run(["g++", "-O3", "-std=c++20", "-march=native", "-mbmi2", str(P / "verify.cpp"), "-o", binary], check=True)
    text = subprocess.check_output([binary], cwd=D, text=True)
a = {}
for line in text.splitlines():
    m = re.match(r"DIM (\d+).*? N (\d+)", line)
    if m and 63 <= int(m[1]) <= 67:
        a[m[1]] = int(m[2])
    if line.startswith("POLAR 65 ADD "):
        a["65"] += int(line.split()[-1])
assert a == json.loads((P / "claims.json").read_text())["lower_bounds"]
print(text, end="")
print(json.dumps({"status": "PASS", "lower_bounds": a}, sort_keys=True))
