"""Exact integer verifier for the complete affine-hyperplane construction."""
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

assert __debug__
P = Path(__file__).resolve().parent
old = P.parent / "dim63-67-signed-overlap"
prev = P.parent / "dim63-68-affine-hyperplane-heads"
witness = {
    "cubic80.txt": (old / "cubic80.txt").read_bytes(),
    "signed_witness.txt": gzip.decompress((old / "signed_witness.txt.gz").read_bytes()),
    "tagged_witness.txt": gzip.decompress((prev / "tagged_witness.txt.gz").read_bytes()),
    "hyperplane_heads.txt": (P / "hyperplane_heads.txt").read_bytes(),
    "joint_witness.txt": (P / "joint_witness.txt").read_bytes(),
}
digest = {
    "cubic80.txt": "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f",
    "signed_witness.txt": "5173fabd124bd23e077d00c956534e7f56543101b28903cbac017ed252c34be7",
    "tagged_witness.txt": "cb51472a7bcbe09f85728b54f74ff1780afaf3a817f4bc442dd7d8b65ba26fb2",
    "hyperplane_heads.txt": "e615223cf61ab75b2a923cc3ecc68b314f6a826c10428b67b5ff10b9da7e21fb",
    "joint_witness.txt": "311ebf0e8de26703fc4f3b790a0f8d81e2f5c183eb154dd2131d530714b409ec",
}
with tempfile.TemporaryDirectory() as d:
    T = Path(d)
    for name, data in witness.items():
        assert hashlib.sha256(data).hexdigest() == digest[name], name
        (T / name).write_bytes(data)
    b = T / "exact"
    subprocess.run(["g++", "-O3", "-march=native", "-std=c++20", "-mbmi2",
                    str(P / "verify.cpp"), "-o", str(b)], check=True)
    out = subprocess.check_output([str(b)], cwd=T, text=True)
    print(out, end="")
claimed = json.loads((P / "claims.json").read_text())
rows = dict(re.findall(r"^DIM (6[5-8]) .* N (\d+)$", out, re.M))
assert claimed == {"status": "PASS", "lower_bounds": {m: int(n) for m, n in rows.items()}}
subprocess.run(["python3", str(P / "local_bounds.py")], check=True)
print(json.dumps(claimed, sort_keys=True))
