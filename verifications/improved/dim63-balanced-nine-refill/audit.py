"""Independent scalar replay of signed packed 63D points."""
from pathlib import Path
import struct,subprocess,tempfile,json,sys,gzip,itertools as it
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from section import setup,D
def audit():
 W,F,H,O=setup()
 data=json.loads((D/'dim63-65-affine-reconstruction/section63.json').read_text())
 P=data['four_support_permutation']
 T=sorted({sum(1<<P[i]for i in (x,y,z,x^y^z))for x,y,z in it.combinations(range(64),3)})
 J=json.loads(gzip.decompress((HERE/'witness63.json.gz').read_bytes()))
 with tempfile.TemporaryDirectory(prefix='63D_scalar_') as p:
  exe=Path(p)/'audit';inp=Path(p)/'witness.bin'
  with inp.open('wb') as out:
   w=lambda *a:out.write(struct.pack('<'+'Q'*len(a),*a))
   w(len(W),len(T),len(O)+len(J),H,0)
   for s in W+T:w(s,0)
   for s,ph in O+J:
    w(s,0,len(ph))
    for v in ph:w(v,0)
  subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'audit.cpp'),'-o',str(exe)],check=True)
  return json.loads(subprocess.check_output([str(exe),str(inp)],text=True))
if __name__=='__main__':print(json.dumps(audit(),sort_keys=True))
