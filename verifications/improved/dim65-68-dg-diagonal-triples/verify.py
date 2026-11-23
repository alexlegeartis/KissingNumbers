#!/usr/bin/env python3
"""Exact atlas refill and three-profile diagonal construction in 65–68D."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run, check_output
from hashlib import sha256
import json,runpy,sys
import verify_diagonal

assert __debug__
here=Path(__file__).resolve().parent
old=here.parent/'dim65-68-dg-retained-atlas'
base65=here.parent/'dim65-dg-three-site-donors'
def pin(dir,name,digest):
    data=(dir/name).read_bytes()
    assert sha256(data).hexdigest()==digest,(dir,name)
    return data
for k,v in {
 'verify.py':'7f2bca666ba844610662f46ef9a1e868f9cafc80385ae59bc423bcfc27890f99',
 'verify_atlas.cpp':'1dac94dd0e8c13244d18c3f97de14786741910d000026334c200c93e607f4fe7',
 'atlas.inc':'95a6b6344cc7a4af5dece3be36e2bb89a311ade757bc9e932334aa9382b5ed85',
 'dg.py':'eefd262c41ae6d93c3407cf41222b0231d40e1e2cf14573b0691c0fece8f0eba',
 'assignment.txt':'12b23c1814c6b4879195c12619b899a77d4a2f614f4d462c79ec65be06bbade7',
 'claims.json':'fd18842c427c0e16ca3be745f8b44f0e2d93583d8ed53a4340ef184444123f26'
}.items():pin(old,k,v)
pin(base65,'verify_gates.py','d2884fbf600218c4753f47a5f65bdd9f98ae32cb5fe23109e998c6b542f9f609')
pin(base65,'triple_assignment.txt','92b4a2dc3284d2ec599c853aa8377d15dc4d85c4ec1a51b4db75f46c3f4a9828')
sys.path.insert(0,str(old))
with redirect_stdout(StringIO()):
    E=runpy.run_path(str(old/'verify.py'))
    runpy.run_path(str(base65/'verify_gates.py'))
assert E['rows']==[(65,40414,180224),(66,45181,360448),(67,54363,540672),(68,68275,2162688)]
assert len(E['C'])==len(set(E['C']))==1024 and len(E['A'])==len(set(E['A']))==128
A=(here/'triple_assignment.txt').read_text().splitlines()
AP=(base65/'triple_assignment.txt').read_text().splitlines()
assert len(A)==23 and A[:-1]==AP[:-1]
heads=set(map(int,A[-1].split()))
assert heads==set(map(int,AP[-1].split()))|{27,115}
marks=[tuple(map(int,line.split()))for line in A[1:-1]]
marked={i for _,*cell in marks for i in cell}
assert len(marked)==63
spare=(set(range(64))-marked).pop()
# The restored ±e_spare do not meet a marked site. They previously cleared
# the retained dense, sparse, four-site and atlas stock. For the new profiles
# the raw dot is at most 1 with squared donor norm >4.
assert spare==7
assert all(n>4*verify_diagonal.D**2 for n in verify_diagonal.N)
phase_delta=json.loads((here/'rephase.json').read_text())
positions=phase_delta['edge_xor_positions']
recovered=phase_delta['recovered_bank_indices']
assert positions==sorted(set(positions)) and len(positions)==68
assert recovered==sorted(set(recovered)) and len(recovered)==51
assert 0<=positions[0]<=positions[-1]<31512
with TemporaryDirectory() as temp:
    T=Path(temp)
    for name,b in E['sources'].items():(T/name).write_bytes(b)
    prev=(T/'assignment.txt').read_text().splitlines()
    prev[0]='46'
    prev[-1]=' '.join(map(str,sorted(heads)))
    (T/'assignment.txt').write_text('\n'.join(prev)+'\n')
    (T/'recover.txt').write_text('\n'.join(map(str,recovered))+'\n')
    w=list((T/'phase.txt').read_text().strip())
    for p in positions:w[p]='1' if w[p]=='0' else '0'
    (T/'phase.txt').write_text(''.join(w)+'\n')
    exe=T/'instrumented'
    run(['g++','-O3','-march=native','-mbmi2','-std=c++20',str(here/'verify_atlas46.cpp'),'-o',str(exe)],check=True)
    receipt=check_output([str(exe)],cwd=T,text=True)
    print(receipt,end='')
    rows=[tuple(map(int,l.split()[1:]))for l in receipt.splitlines()if l.startswith('ATLAS ')]
    assert rows==[(65,40414,188416),(66,45181,376832),(67,54363,565248),(68,68275,2260992)]
    crt=T/'cross'
    run(['g++','-O3','-std=c++20',str(here/'cross.cpp'),'-o',str(crt)],check=True)
    run([str(crt),str(T/'cubic80.txt'),str(here/'triple_assignment.txt'),str(T/'hyperplane_heads.txt')],check=True)

    mon=lambda u:sum(1<<i for i in range(64)if i&u==u)
    quads=[mon(u)for u in range(64)if u.bit_count()<=2]
    cub=[mon(u)for u in range(64)if u.bit_count()==3]
    def vx(v):
        w=0
        for i,h in enumerate(cub):
            if v>>i&1:w^=h
        return w
    cm=[vx(int(d))for d in (T/'cubic80.txt').read_text().split()]
    assert len(cm)==len(set(cm))==80
    triples=[(cm[j],sum(1<<i for i in cell))for j,*cell in marks]
    drop=[];nt=0;odd=0
    for l in (T/'tag3.txt').read_text().splitlines():
        i,b,sz,z,*V=map(int,l.split())
        assert b.bit_count()==14 and len(V)==sz and V[:2]==[0,65535]
        assert all(not ((u&b).bit_count()&1) for u in quads)
        sp=(z>>1)&1
        nt+=1;odd+=sp
        # With two marked donor tail entries, both signs are equal.
        # If any exact magnitude gate fails, a marked large negative
        # site lies in this bank. A mismatch follows from the total
        # parity unless dot(cubic,b) differs from the sparse parity.
        if any((b&m).bit_count()>=2 and ((c&b).bit_count()&1)!=sp for c,m in triples):
            assert sp
            drop.append(b)
    assert nt==4767 and odd==3 and not drop
    print('PASS all 4,767 two-tail banks, no deletions')
oldn=E['out']['lower_bounds']
bank=len(E['A'])*len(E['C'])
gain=10*bank
ans={str(d):oldn[str(d)]+gain+2+(heads1-heads0)
     for (d,_,heads1),(old_d,_,heads0) in zip(rows,E['rows']) if d==old_d}
assert ans=={'65':424177492,'66':434139558,'67':453138696,'68':483338010}
out={'status':'PASS','lower_bounds':ans}
assert out==json.loads((here/'claims.json').read_text())
print(json.dumps(out,sort_keys=True))
