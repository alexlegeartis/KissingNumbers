#!/usr/bin/env python3
"""Independently reconstruct and certify the self-sparse 62–66D family."""
from itertools import combinations, product
from pathlib import Path
import hashlib, json, os, subprocess, tempfile
assert __debug__  # Exact replay must not suppress assertions.
H=Path(__file__).resolve().parent
p=H/'witness.txt'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f'
C=list(map(int,p.read_text().split()))
assert len(C)==len(set(C))==80

with tempfile.TemporaryDirectory() as tmp:
    exe={}
    for name in ('gray','plotkin','verify_independent','tails'):
        exe[name]=str(Path(tmp)/name)
        subprocess.run(['g++','-std=c++17','-O3','-march=native',str(H/f'{name}.cpp'),'-o',exe[name]],check=True)
    # Independently check the true length-64 distance with a second full enumerator.
    r=subprocess.check_output([exe['verify_independent'],str(p)],cwd=tmp,text=True)
    assert 'PASS independent' in r
    print(r.strip())
    # Distinct Gray and Plotkin routes independently reconstruct the entire sparse stock.
    gray=subprocess.check_output([exe['gray'],str(p),'supports64.bin'],cwd=tmp,text=True)
    plotkin=subprocess.check_output([exe['plotkin'],str(p)],cwd=tmp,text=True)
    print(gray.strip())
    print(plotkin.strip())
    a=(Path(tmp)/'supports64.bin').read_bytes()
    U=[int.from_bytes(a[j:j+8],byteorder=os.sys.byteorder)for j in range(0,len(a),8)]
    V=list(map(int,(Path(tmp)/'supports64.txt').read_text().split()))
    assert U==V==sorted(set(V)) and len(V)==37996
    tails=subprocess.check_output([exe['tails'],str(p)],cwd=tmp,text=True)
    print(tails.strip())
    assert 'tail support counts 37996 3904 3904 3904' in tails
assert all(w.bit_count()==16 for w in V)

# All affine planes in F2^6; any distinct two have at most 2 common sites.
P={sum(1<<x for x in (a,b,c,a^b^c)) for a,b,c in combinations(range(64),3)}
assert len(P)==10416 and all(x.bit_count()==4 for x in P)
# Explicit local [16,11,4] parity check. Its dual has no word of weight 2.
checks=[65535]+[sum(1<<i for i in range(16)if (i>>j)&1)for j in range(4)]
D={0}
for h in checks:D|={h^s for s in D}
assert len(D)==32 and {s.bit_count()for s in D}=={0,8,16}
signs=[s for s in range(65536)if all((s&h).bit_count()%2==0 for h in checks)]
assert len(signs)==2048 and min(s.bit_count() for s in signs if s)==4

# 63D equality is an exact linear section of the entire 64D code.
a,b=0,20;m=(1<<a)|(1<<b)
bitrows=[sum((int(x&s==s)<<i)for i,x in enumerate((a,b)))for s in range(64)if s.bit_count()<=2]
R={0}
for z in bitrows:R|={z^x for x in R}
assert len(R)==4
s0=sum(not m&w for w in V);s2=sum(m&w==m for w in V)
q0=sum(not m&w for w in P);q2=sum(m&w==m for w in P)
assert (s0,s2,q0,q2)==(21361,2363,9145,31)
n63=(80<<21)+2048*s0+1024*s2+16*q0+8*q2+124

# 62D: substitute hexagon labels for four marked binary coordinates,
# including a fixed zero in the fourth; all ranks checked on degree <=2.
marks=(0,20,41,61);m=sum(1<<x for x in marks)
R={0}
for s in range(64):
    if s.bit_count()>2:continue
    z=sum(int(x&s==s)<<j for j,x in enumerate(marks));R|={z^x for x in R}
assert len(R)==16
labels=(0,1,3,7,6,4)
heads=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
for (i,(x,y)),(j,(u,v)) in product(zip(labels,heads),repeat=2):
    assert 60+x*x+3*y*y==64
    assert 60-2*(16-(i^j).bit_count())+x*u+3*y*v <= 32
    # Cross with amplitude-4 affine planes containing neither, one or both head sites.
    for a,b in product((0,1),repeat=2):
        k=32-4*(4-a-b+a*abs(x))
        assert k>=0 and k*k>=48*b*y*y
    assert 8*abs(x)<=32 and 192*y*y<=32*32
# We keep two of the marked coordinate names for the two physical head axes.
# The plane supports live on the 62 chosen names; the sparse vectors use only the 60 body names.
removed=(1<<41)|(1<<61)
plane62=sum(not w&removed for w in P)
s62=sum(not w&m for w in V)
assert (plane62,s62)==(9145,12005)
n62=6*(80<<18)+s62*2048+plane62*16+124
n64=(80<<22)+len(V)*2048+len(P)*16+128
actual={'status':'PASS','supports64':len(V),'supports62':s62,'lower_bounds':
        {str(d):n for d,n in [(62,n62),(63,n63),(64,n64),(65,n64+2),(66,n64+4+3904*2048),(67,n64+6+3*3904*2048)]}}
assert actual==json.loads((H/'claims.json').read_text())
print(json.dumps(actual,indent=2))
