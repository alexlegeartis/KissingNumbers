#!/usr/bin/env python3
"""Exact replay of the 90D constant-rail marked QR16 witness."""
from pathlib import Path
import base64, hashlib, os, struct, subprocess, tempfile, zlib
here=Path(__file__).resolve().parent
phase=(here/"phase.b64").read_bytes()
assert hashlib.sha256(phase).hexdigest()=="80d6285fbb9bebe81d4014b7dbce7df0393e495537b82c298ba6971bda464091"
values=[x for x, in struct.iter_unpack("<I",zlib.decompress(base64.b64decode(phase)))]
assert len(values)==len(set(values))==7906 and all(0<=n<3**15 for n in values)
def be3(n):
    s=""
    for _ in range(15): s=str(n%3)+s; n//=3
    assert not n
    return s
K=("\n".join(be3(x) for x in values)+"\n").encode()
assert hashlib.sha256(K).hexdigest()=="220a8d800eaaf391898b424e5e27f54fe957c4f0548d6ce980e9b488637ab369"
w=base64.b64decode((here/"witness.b64").read_bytes())
assert hashlib.sha256(w).hexdigest()=="c5b0e90a8554d0a14a9cb6191eb795d359f6d83d5928feea6950972879320186"
with tempfile.TemporaryDirectory() as d:
    d=Path(d); p=d/"witness.bin"; p.write_bytes(w)
    tool=d/"audit"; subprocess.run([os.environ.get("CXX","g++"),"-std=c++20","-O3",str(here/"verify.cpp"),"-o",str(tool)],check=True)
    completed=subprocess.run([str(tool),str(p)],input=K,stdout=subprocess.PIPE,check=True)
    line=completed.stdout.decode().strip()
    assert line.startswith("AUDIT PASS ")
    print(line)
