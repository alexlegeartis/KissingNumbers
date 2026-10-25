"""Verify the 73-dimensional lattice and a 6,218,188,144-point lower bound.

NumPy and Python's standard library suffice. Use --rebuild with g++ installed
to rerun the three complete interval enumerations (about one minute serially).
The parent Gamma72 minimum and kissing number use Nebe's published theorem.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import numpy as np

HERE=Path(__file__).resolve().parent
MASK=(1<<24)-1


def rational_ldl(g):
    n=len(g)
    mu=[[F(i==j) for j in range(n)] for i in range(n)]
    d=[]
    for i in range(n):
        d.append(F(int(g[i,i]))-sum(mu[i][k]**2*d[k] for k in range(i)))
        assert d[-1]>0
        for j in range(i+1,n):
            mu[j][i]=(int(g[j,i])-sum(mu[j][k]*mu[i][k]*d[k] for k in range(i)))/d[i]
    return d,mu


def bits(row):
    return sum((int(v)&1)<<i for i,v in enumerate(row))


def bit_table(matrix):
    rows=[bits(row) for row in matrix]
    result=np.zeros((3,256),dtype=np.uint32)
    for part in range(3):
        for x in range(256):
            for j in range(8):
                if x>>j&1: result[part,x]^=rows[8*part+j]
    return result


def apply(table,x):
    return table[0,x&255]^table[1,(x>>8)&255]^table[2,x>>16]


def count_triples(lists,a,b,h):
    """Independent count using full 24-bit labels, without search score files."""
    ta,tb=bit_table(a),bit_table(b)
    groups=[]
    for raw in lists[:2]:
        z=raw&MASK
        keys=((raw>>24).astype(np.uint64)<<24)|apply(tb,z)
        order=np.argsort(keys);keys,z=keys[order],z[order]
        cut=np.r_[0,np.flatnonzero(keys[1:]!=keys[:-1])+1,len(keys)]
        grouped={}
        for i,j in zip(cut[:-1],cut[1:]):
            key=int(keys[i]);grouped[(key>>24,key&MASK)]=z[i:j]
        groups.append(grouped)
    third=np.full(1<<24,255,dtype=np.uint8)
    third[lists[2]&MASK]=(lists[2]>>24).astype(np.uint8)
    shorter,nearest,pairs=0,0,0
    for (c1,label),z1 in groups[0].items():
        az1=apply(ta,z1)
        for c2 in range(2,6):
            if c1+c2>9:continue
            z2=groups[1].get((c2,label))
            if z2 is None:continue
            keys=(az1[:,None]^z2[None,:]^h).ravel()
            cost=third[keys].astype(np.int32)+c1+c2
            shorter+=int(np.count_nonzero(cost<=9))
            nearest+=int(np.count_nonzero(cost==11))
            pairs+=len(keys)
    assert shorter==0
    assert nearest==6008
    return {'shorter_triples':shorter,'nearest_with_all_block_costs_at_most_5':nearest,'checked_pairs':pairs}


def rebuild(lists,qs):
    compiler=shutil.which('g++')
    if not compiler:raise RuntimeError('g++ is required for --rebuild')
    with tempfile.TemporaryDirectory(prefix='lattice73-') as tmp:
        tmp=Path(tmp)
        exe=tmp/'enumerate.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'enumerate_leech.cpp'),'-o',str(exe)],check=True)
        for i in range(3):
            q=tmp/f'q{i}.txt';q.write_text(' '.join(map(str,qs[i]))+'\n')
            out=tmp/f'list{i}'
            subprocess.run([str(exe),str(HERE/'leech_ldl.txt'),str(q),'20','3600',str(out),'interval'],check=True)
            metadata=json.loads(out.with_suffix('.json').read_text())
            assert metadata['complete'] and metadata['certified_intervals']
            assert np.array_equal(np.fromfile(out.with_suffix('.bin'),dtype=np.uint32),lists[i])
        target_exe=tmp/'targets.exe'
        subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'enumerate_targets.cpp'),'-o',str(target_exe)],check=True)
        out=tmp/'completion'
        subprocess.run([str(target_exe),str(HERE/'leech_ldl.txt'),str(HERE/'completion.input.txt'),str(out)],check=True)
        assert out.with_suffix('.vectors.bin').read_bytes()==(HERE/'completion.vectors.bin').read_bytes()


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--rebuild',action='store_true')
    args=parser.parse_args()
    started=time.time()
    data=np.load(HERE/'certificate.npz')
    f,a,b,t,g,x,u=[data[k].astype(np.int64) for k in ('gram','a','b','t','gamma_gram','heads','centre')]
    published=np.array(json.loads((HERE/'catalogue_gram72.json').read_text()),dtype=np.int64)
    assert np.array_equal(g,published)
    eye=np.eye(24,dtype=np.int64)
    assert np.array_equal(a+b,eye)
    assert np.array_equal(a@b,2*eye)
    assert np.array_equal(a@f@a.T,2*f)
    z=np.zeros((24,24),dtype=np.int64)
    assert np.array_equal(t,np.block([[eye,eye,b],[z,a,a],[z,z,2*eye]]))
    assert np.array_equal(t@np.kron(np.eye(3,dtype=np.int64),f)@t.T,2*g)
    assert np.array_equal(g,g.T) and not (g.diagonal()%2).any()
    d,mu=rational_ldl(f)
    det=F(1)
    for q in d:det*=q
    assert det==1 and not (f.diagonal()%2).any()
    # det(a)^2=2^24; det(t)^2=2^72, so det(g)=1.
    vals=(HERE/'leech_ldl.txt').read_text().split()
    assert int(vals[0])==24
    assert np.array_equal(np.array(list(map(int,vals[1:577]))).reshape(24,24),f)
    interval=list(map(float,vals[577:]));exact=d+[v for row in mu for v in row]
    assert len(interval)==2*len(exact)
    for j,q in enumerate(exact):
        assert F.from_float(interval[2*j])<=q<=F.from_float(interval[2*j+1])
    qs=[a[0],a[1],a[0]+a[1]]
    lists=[data[f'list{i}'] for i in range(3)]
    for raw in lists:
        assert len(raw)==1231792
        assert dict(zip(*np.unique(raw>>24,return_counts=True)))=={2:48,3:4096,4:97152,5:1130496}
        assert len(np.unique(raw&MASK))==len(raw)
    if args.rebuild:rebuild(lists,qs)

    assert x.shape==(6272,72) and np.array_equal(x[0],u)
    assert len(np.unique(x,axis=0))==6272
    assert not ((x-u)%2).any()
    assert np.all(np.einsum('ij,jk,ik->i',x,g,x)==22)
    keys={tuple(row) for row in x}
    assert all(tuple(-row) in keys for row in x)
    assert u@g@u==22
    # Derive the label h from the given centre, rather than hardcoding it.
    w=(u@t).reshape(3,24)
    assert not ((w-np.array(qs))%2).any()
    labels=np.array([bits(row) for row in (w-np.array(qs))//2],dtype=np.uint32)
    ta,tb=bit_table(a),bit_table(b)
    assert apply(tb,labels[0])==apply(tb,labels[1])
    h=int(labels[2]^apply(ta,labels[0])^labels[1])
    assert h==12625758
    counts=count_triples(lists,a,b,h)

    # Cross-check the source of the 264 new representatives in exact arithmetic.
    wx=(x@t).reshape(-1,3,24)
    block_norms=np.einsum('bki,ij,bkj->bk',wx,f,wx)
    assert not (block_norms%4).any()
    costs=block_norms//4
    assert np.all(costs.sum(axis=1)==11) and costs.min()==2
    assert int(np.count_nonzero(costs.max(axis=1)<=5))==6008
    assert int(np.count_nonzero(costs.max(axis=1)>5))==264
    patterns=Counter(tuple(sorted(map(int,row))) for row in costs)

    h73=np.array(json.loads((HERE/'gram73_twice.json').read_text()),dtype=np.int64)
    expected=np.zeros((73,73),dtype=np.int64)
    expected[:72,:72]=2*g
    expected[:72,72]=g@u;expected[72,:72]=u@g;expected[72,72]=16
    assert np.array_equal(h73,expected)
    # Schur complement of G73 is 8-22/4=5/2: positive definite, det 5/2.
    schur=F(8)-F(int(u@g@u),4)
    assert schur==F(5,2)
    # Every layer |k|>=2 has squared height >=10; layers +/-1 use u/2+Gamma.
    assert 4*schur>8 and F(22,4)+schur==8
    added=np.vstack([np.column_stack(((x-u)//2,np.ones(len(x),dtype=np.int64))),
                     np.column_stack(((x+u)//2,-np.ones(len(x),dtype=np.int64)))])
    assert len(np.unique(added,axis=0))==12544
    assert np.all(np.einsum('ij,jk,ik->i',added,h73,added)==16)
    result={'status':'verified','dimension':73,'minimum_squared_norm':8,
            'determinant':'5/2','certified_coset_vectors':len(x),
            'kissing_number_lower_bound':6218175600+len(added),
            'previous_repository_bound':6218187616,'gain':528,
            'block_patterns':{str(k):v for k,v in sorted(patterns.items())},
            'minimum_check':counts,'enumerations_rebuilt':args.rebuild,
            'seconds':time.time()-started,
            'scope':'Nebe parent minimum and shell count are literature inputs; no global priority assertion.'}
    assert result['kissing_number_lower_bound']==6218188144
    report_name='verification_rebuild.json' if args.rebuild else 'verification.json'
    (HERE/report_name).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
