#!/usr/bin/env python3
"""Standalone exact 90-dimensional QR18 certificate replay."""
from pathlib import Path
import base64, gzip, hashlib, os, struct, subprocess, tempfile, zlib
p=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
x=(p/"witness.bin.gz").read_bytes()
assert sha(x)=="b8b3c0060a5b7b51f4dd1aeef5e735e8dde6d25c3d375043d7027b4f43070209"
w=gzip.decompress(x)
assert sha(w)=="ba1695bee7fa7a3a3c3fc5b600e580b65f720579e17383d4545167600ed18c41"
a=(p/"phase.b64").read_bytes()
assert sha(a)=="80d6285fbb9bebe81d4014b7dbce7df0393e495537b82c298ba6971bda464091"
v=[n for n, in struct.iter_unpack("<I",zlib.decompress(base64.b64decode(a)))]
assert len(v)==len(set(v))==7906 and all(0<=n<3**15 for n in v)
def fmt(n):
 s=""
 for j in range(15):
  s=str(n%3)+s;n//=3
 assert not n
 return s
ph=("\n".join(map(fmt,v))+"\n").encode()
assert sha(ph)=="220a8d800eaaf391898b424e5e27f54fe957c4f0548d6ce980e9b488637ab369"
with tempfile.TemporaryDirectory() as t:
 t=Path(t);(t/"witness.bin").write_bytes(w)
 command=[os.environ.get("CXX","g++"),"-std=c++20","-O3",str(p/"verify.cpp"),"-o",str(t/"verify")]
 subprocess.run(command,check=True)
 r=subprocess.run([str(t/"verify"),str(t/"witness.bin")],input=ph,check=True,stdout=subprocess.PIPE,text=False)
 assert r.stdout.decode().strip()=="AUDIT PASS count=758 old_shifts=82452 near=269246 exact_total=21956420880"
 print(r.stdout.decode().strip())

