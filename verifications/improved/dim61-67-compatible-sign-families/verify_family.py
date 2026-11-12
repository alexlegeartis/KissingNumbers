"""Reconstruct finite sign families and compare two exact arithmetic paths."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,itertools,json,subprocess,tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'verifications/improved/dim61-68-95-96-code-layers/source'

def rref(rows,n):
    rows=list(rows);k=0;piv=[]
    for i in range(n):
        j=next((j for j in range(k,len(rows)) if rows[j]>>i&1),None)
        if j is None:continue
        rows[k],rows[j]=rows[j],rows[k]
        for j in range(len(rows)):
            if j!=k and rows[j]>>i&1:rows[j]^=rows[k]
        piv.append(i);k+=1
    assert not any(rows[k:])
    return rows[:k],piv

def generator():
    exponents=(37,36,34,33,32,27,25,24,22,21,19,18,15,11,10,8,7,5,3,0)
    g=sum(1<<i for i in exponents)
    z=(1<<73)|1
    while z.bit_length()>=g.bit_length():z^=g<<(z.bit_length()-g.bit_length())
    assert z==0
    rows=[g<<i for i in range(36)];deleted=(21,37,38,49,57,63,67,70,72)
    for i in deleted:
        j=next((j for j,x in enumerate(rows) if x>>i&1),None)
        if j is None:continue
        p=rows.pop(j);rows=[x^p if x>>i&1 else x for x in rows]
    assert len(rows)==28
    keep=[i for i in range(73) if i!=0 and i not in deleted]
    rows=[sum(((x>>c)&1)<<j for j,c in enumerate(keep)) for x in rows]
    rows=[x|((x.bit_count()%2)<<63) for x in rows]
    rows=rref(rows,64)[0];assert len(rows)==28
    assert all((a&b).bit_count()%2==0 for a in rows for b in rows)
    return rows

def span(rows):
    ans=[0]
    for row in rows:ans += [x^row for x in ans]
    return ans

def expand(s,w):
    return sum(((w>>j)&1)<<c for j,c in enumerate(i for i in range(64) if s>>i&1))

def code(row,odd):
    s=row['support'];h=row['checks'];t=row['shift']
    assert isinstance(s,int) and 0<s<1<<64 and s.bit_count()==16
    assert len(h)==5 and h[0]==65535 and all(isinstance(x,int) and 0<=x<65536 for x in h)
    assert isinstance(t,int) and 0<=t<65536 and (not odd or t.bit_count()%2==1)
    cols=[sum(((a>>i)&1)<<j for j,a in enumerate(h)) for i in range(16)]
    assert sorted(cols)==list(range(1,32,2))
    dual=span(h);assert len(set(dual))==32 and {v.bit_count() for v in dual}=={0,8,16}
    return dict(s=s,dual={expand(s,v) for v in dual},phase=expand(s,t))

def popcount(x):
    x=x-((x>>np.uint64(1))&np.uint64(0x5555555555555555))
    x=(x&np.uint64(0x3333333333333333))+((x>>np.uint64(2))&np.uint64(0x3333333333333333))
    x=(x+(x>>np.uint64(4)))&np.uint64(0x0f0f0f0f0f0f0f0f)
    return ((x*np.uint64(0x0101010101010101))>>np.uint64(56)).astype(np.uint8)

def check_pair(a,b):
    overlap=a['s']&b['s'];assert overlap.bit_count()==10
    common=[h for h in a['dual']&b['dual'] if h.bit_count()==8 and h&overlap==h]
    assert len(common)==1
    h=common[0];assert (a['phase']&h).bit_count()%2!=(b['phase']&h).bit_count()%2

def check(path,cpp=True):
    if not __debug__:raise RuntimeError('Run without -O.')
    raw=path.read_bytes();cert=json.loads(raw)
    archive=SOURCE/'dim61_96/code64_supports.npz'
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==cert['source_archive_sha256']
    with np.load(archive,allow_pickle=False) as z:
        old=list(map(int,z['words']));saved=[sum(int(v)<<i for i,v in enumerate(r)) for r in z['generator']]
    gen=generator();assert rref(saved,64)[0]==gen
    assert len(old)==len(set(old))==30828 and all(s.bit_count()==16 for s in old)
    deleted=cert['removed_coordinates'];assert len(deleted)==len(set(deleted)) and all(0<=x<64 for x in deleted)
    removed=sum(1<<i for i in deleted);odd=cert['odd_total_parity'];assert isinstance(odd,bool)
    if removed:assert deleted==[24,42,53,55] and odd
    kept=[s for s in old if not s&removed]
    assert len(kept)==cert['original_supports_retained']
    new=[code(row,odd) for row in cert['new']]
    modified={row['support']:code(row,odd) for row in cert['old']}
    assert len(modified)==len(cert['old']) and set(modified)<=set(kept)
    assert len({a['s'] for a in new})==len(new) and not {a['s'] for a in new}&set(old)
    assert all(not a['s']&removed and not((a['s']&g).bit_count()%2) for a in new for g in gen)
    oa=np.array(kept,dtype=np.uint64);na=np.array([a['s'] for a in new],dtype=np.uint64)
    oc=np.array([[(s>>i)&1 for i in range(64)] for s in kept],dtype=np.float32)
    nc=np.array([[(a['s']>>i)&1 for i in range(64)] for a in new],dtype=np.float32)
    repairs_old=repairs_new=0;touched=set()
    for i,a in enumerate(new):
        packed=popcount(oa&np.uint64(a['s']));physical=oc@nc[i]
        assert np.array_equal(packed,physical) and int(packed.max())<=10
        for j in np.flatnonzero(packed>8):
            b=kept[int(j)];assert b in modified;check_pair(a,modified[b]);touched.add(b);repairs_old+=1
        if i:
            packed=popcount(na[:i]&np.uint64(a['s']));physical=nc[:i]@nc[i]
            assert np.array_equal(packed,physical) and int(packed.max())<=10
            for j in np.flatnonzero(packed>8):check_pair(a,new[int(j)]);repairs_new+=1
    assert touched==set(modified)
    scalar=None
    if cpp:
        with tempfile.TemporaryDirectory(prefix='parity-family-') as td:
            td=Path(td);exe=td/'check';text=td/'input.txt'
            rows=[f'{len(old)} {len(new)} {len(modified)} {len(gen)} {int(odd)} {removed}',
                  ' '.join(map(str,gen)),' '.join(map(str,old))]
            for group in (cert['old'],cert['new']):
                rows.extend(' '.join(map(str,[x['support'],*x['checks'],x['shift']])) for x in group)
            text.write_text('\n'.join(rows)+'\n')
            subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'check_family.cpp'),'-o',str(exe)],check=True)
            result=subprocess.run([str(exe),str(text)],check=True,capture_output=True,text=True)
            scalar=json.loads(result.stdout)
            assert scalar['old_new_repairs']==repairs_old and scalar['new_new_repairs']==repairs_new
            assert scalar['old_new_pairs']==len(kept)*len(new) and scalar['new_new_pairs']==len(new)*(len(new)-1)//2
    return dict(status='PASS',certificate_sha256=hashlib.sha256(raw).hexdigest(),
                retained_old_supports=len(kept),new_supports=len(new),modified_old=len(modified),
                old_new_pairs=len(kept)*len(new),new_new_pairs=len(new)*(len(new)-1)//2,
                old_new_repairs=repairs_old,new_new_repairs=repairs_new,scalar_cpp=scalar)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('certificate',type=Path);parser.add_argument('--no-cpp',action='store_true');args=parser.parse_args()
    print(json.dumps(check(args.certificate,not args.no_cpp),indent=2))
