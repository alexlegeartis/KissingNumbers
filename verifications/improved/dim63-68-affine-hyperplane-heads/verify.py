import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

assert __debug__
P = Path(__file__).resolve().parent
PRE = P.parent / "dim63-67-signed-overlap"
FILES = {
    "cubic80.txt": "5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f",
    "signed_witness.txt": "5173fabd124bd23e077d00c956534e7f56543101b28903cbac017ed252c34be7",
    "tagged_witness.txt": "cb51472a7bcbe09f85728b54f74ff1780afaf3a817f4bc442dd7d8b65ba26fb2",
    "hyperplane_heads.txt": "d446cf0b00ea3769d08165921775fb74ad7761dfc460913a7b9561e75f1a63cc",
    "joint_witness.txt": "ff1375eff3010defcdd57145c402d17664adac39eea8d1fcddc084e08d67b87e",
    "tall.txt": "734bb49dadfeffbaa0404cae35d4f904a916f0f56be9f7a3e82e1cb5a0201c94",
}
PAYLOADS = {
    "cubic80.txt": (PRE / "cubic80.txt").read_bytes(),
    "signed_witness.txt": gzip.decompress((PRE / "signed_witness.txt.gz").read_bytes()),
    "tagged_witness.txt": gzip.decompress((P / "tagged_witness.txt.gz").read_bytes()),
    "hyperplane_heads.txt": (P / "hyperplane_heads.txt").read_bytes(),
    "joint_witness.txt": (P / "joint_witness.txt").read_bytes(),
    "tall.txt": (P / "tall.txt").read_bytes(),
}
counts = {}
with tempfile.TemporaryDirectory() as d:
    T = Path(d)
    for name, payload in PAYLOADS.items():
        assert hashlib.sha256(payload).hexdigest() == FILES[name], name
        (T / name).write_bytes(payload)
    def run(source):
        b = T / Path(source).stem
        subprocess.run(["g++", "-O3", "-std=c++20", "-march=native", "-mbmi2",
                        str(P / source), "-o", str(b)], check=True)
        return subprocess.check_output([str(b)], cwd=T, text=True)
    a = run("verify.cpp")
    b = run("polar.cpp")
    m = re.search(r"^DIM 65 .* N (\d+)$", a, re.M)
    assert m and int(m[1]) == 418478978
    p = re.search(r"^65 (\d+)$", b, re.M)
    assert p
    counts["65"] = int(p[1]) + int(m[1]) - 418347906
    for n in (63, 64, 66, 67, 68):
        m = re.search(rf"^DIM {n} (?:.* )?N (\d+)(?: at .*|)$", a, re.M)
        assert m
        counts[str(n)] = int(m[1])
    print(a, end="")
    print(b.splitlines()[1])
claims = json.loads((P / "claims.json").read_text())
assert claims == {"status": "PASS", "lower_bounds": counts}
print(json.dumps(claims, sort_keys=True))
