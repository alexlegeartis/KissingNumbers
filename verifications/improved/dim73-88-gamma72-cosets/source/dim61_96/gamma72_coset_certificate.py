#!/usr/bin/env python3
"""Independent verification of the Gamma72 nearest-coset constructions, 73--77.

Run from the repository root:
    python research/dim61_96/gamma72_coset_certificate.py

The 24-dimensional lists were exhaustively enumerated using a rational LDL
decomposition enclosed by outward rounded binary64 intervals.  This verifier
checks the enclosures and completed-search metadata, independently checks the
absence of shorter triples at every coset actually used, and reconstructs and
checks explicit norm-22 representatives from the selected-vector exports.

No assertion of completeness is made for the retained norm-22 shell: these
are certified subsets of nearest vectors.  The full Gamma72 minimum 8 and
its 6,218,175,600-vector shell are from Nebe's extremality theorem.

Geometric lemma.  If S is any subset of nearest vectors of c+L, all of norm
r<=6, then (s,+/-sqrt(8-r)) coexist with the norm-8 shell of L.  For distinct
cosets, squared coset distance d and tail distance e suffice when d+e>=8.

Sources:
  https://arxiv.org/abs/1008.2862
  https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb72.html
"""
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
import time

import numpy as np
import sympy as sp

from gamma72_cosets import HERE, SOURCE, rational_ldl

BASE = 6218175600
HEADS = {12625758: 6008, 4915087: 6000, 6896114: 5998}
PREVIOUS = {73: 6218179814, 74: 6218192454, 75: 6218209308,
            76: 6218243016, 77: 6218285152}
MASK = (1 << 24)-1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(v):
    return sum((int(x)&1)<<i for i,x in enumerate(v))


def unbits(x):
    return np.array([(int(x)>>i)&1 for i in range(24)], dtype=np.int64)


def table(matrix):
    rows = [bits(row) for row in matrix]
    result = np.zeros((3,256),dtype=np.uint32)
    for part in range(3):
        for value in range(256):
            acc = 0
            for j in range(8):
                if (value>>j)&1:
                    acc ^= rows[8*part+j]
            result[part,value] = acc
    return result


def apply(t, x):
    return t[0,x&255] ^ t[1,(x>>8)&255] ^ t[2,x>>16]


def check_intervals(f):
    values = (HERE/'gamma72_leech_ldl.txt').read_text().split()
    assert int(values[0]) == 24
    assert np.array_equal(np.array(list(map(int,values[1:577]))).reshape(24,24), f)
    doubles = list(map(float,values[577:]))
    d, mu = rational_ldl(f)
    exact = d + [x for row in mu for x in row]
    assert len(doubles) == 2*len(exact)
    for k,x in enumerate(exact):
        assert Fraction.from_float(doubles[2*k]) <= x <= Fraction.from_float(doubles[2*k+1])


def load_lists(power, hashes, qfamily):
    tag='0_1'+('_d%d'%power if power else '')
    arrays=[]
    for block in range(3):
        input_path=HERE/('gamma72_leech_%s_%d.input.txt'%(tag,block))
        actual_q=np.loadtxt(input_path,dtype=np.int64)
        assert np.array_equal(actual_q,qfamily[block])
        hashes[input_path.name]=sha(input_path)
        prefix=HERE/('gamma72_leech_%s_%d_r5_cert'%(tag,block))
        metadata=json.loads(prefix.with_suffix('.json').read_text())
        assert metadata['complete'] and metadata['certified_intervals']
        assert metadata['radius_squared_leech'] == 5
        path=prefix.with_suffix('.bin')
        records=np.fromfile(path,dtype=np.uint32)
        assert len(records) == 1231792
        assert dict(zip(*np.unique(records>>24,return_counts=True))) == {
            2:48,3:4096,4:97152,5:1130496}
        assert len(np.unique(records&MASK)) == len(records)
        hashes[path.name]=sha(path)
        arrays.append(records)
    return arrays


def no_shorter_triples(records, a, target_h):
    """Independent direct lookup, with full 24-bit N labels (no compact labels).

    Each block has minimum cost 2.  A total cost <=9 therefore uses only
    block costs <=5 and pair costs <=7.  These are exhaustively checked.
    """
    ta,tb=table(a),table(np.eye(24,dtype=np.int64)-a)
    groups=[]
    for raw in records[:2]:
        z=raw&MASK
        keys=((raw>>24).astype(np.uint64)<<24) | apply(tb,z)
        order=np.argsort(keys)
        keys,z=keys[order],z[order]
        az=apply(ta,z)
        starts=np.r_[0,np.flatnonzero(keys[1:]!=keys[:-1])+1,len(keys)]
        grouped={}
        for i,j in zip(starts[:-1],starts[1:]):
            key=int(keys[i]);grouped[(key>>24,key&MASK)]=(z[i:j],az[i:j])
        groups.append(grouped)
    third=np.full(1<<24,255,dtype=np.uint8)
    third[records[2]&MASK]=(records[2]>>24).astype(np.uint8)
    found={int(h):0 for h in target_h}
    pairs=0
    for (c1,nlabel),(z1,az1) in groups[0].items():
        for c2 in range(2,min(5,7-c1)+1):
            second=groups[1].get((c2,nlabel))
            if second is None:
                continue
            pairkeys=(az1[:,None]^second[0][None,:]).reshape(-1)
            pairs+=len(pairkeys)
            for h in found:
                found[h]+=int(np.count_nonzero(third[pairkeys^h] <= 9-c1-c2))
    assert all(count==0 for count in found.values()),found
    return {'checked_pair_products':pairs,'shorter_triple_counts':found}


def to_catalogue(w, a, b):
    r1=w[:,0,:]
    twice_r2=(w[:,1,:]-r1)@b
    twice_r3=w[:,2,:]-w[:,1,:]+r1@a
    assert not (twice_r2%2).any() and not (twice_r3%2).any()
    return np.hstack([r1,twice_r2//2,twice_r3//2])


def check_heads(a,b,f,g72,hashes):
    selected=[]
    source=[]
    for block in range(3):
        p=HERE/('gamma72_export_block%d.selected.bin'%block)
        rows=np.fromfile(p,dtype=np.int32).reshape(-1,25).astype(np.int64)
        assert (np.diff(rows[:,0])>0).all()
        selected.append(rows)
        cert=HERE/('gamma72_leech_0_1_%d_r5_cert.bin'%block)
        assert sha(HERE/('gamma72_export_block%d.bin'%block)) == sha(cert)
        source.append(np.fromfile(cert,dtype=np.uint32))
        hashes[p.name]=sha(p)
    result={}
    stored={}
    for h,count in HEADS.items():
        triples=np.fromfile(HERE/('gamma72_head_%d.triples.bin'%h),dtype=np.uint32).reshape(-1,3)
        assert len(triples)==count
        blocks=[]
        for block in range(3):
            rows=selected[block]
            positions=np.searchsorted(rows[:,0],triples[:,block])
            assert np.array_equal(rows[positions,0],triples[:,block])
            w=rows[positions,1:]
            q=a[0] if block==0 else a[1] if block==1 else a[0]+a[1]
            assert not ((w-q)%2).any()
            norms=np.einsum('ij,jk,ik->i',w,f,w)
            raw=source[block][triples[:,block]]
            assert np.array_equal(norms,4*(raw>>24))
            z=(w-q)//2
            assert np.array_equal(np.array([bits(x) for x in z],dtype=np.uint32),raw&MASK)
            blocks.append(w)
        w=np.stack(blocks,axis=1)
        reps=to_catalogue(w,a,b)
        norms=np.einsum('ij,jk,ik->i',reps,g72,reps)
        assert (norms==22).all()
        assert len(np.unique(reps,axis=0))==count
        assert not ((reps-reps[0])%2).any()
        keys={row.tobytes() for row in reps}
        assert all((-row).tobytes() in keys for row in reps)
        # Exact arithmetic in binary64: all terms and all partial sums are
        # integers below 2^53.  The bound checks this before BLAS is used.
        rg=reps@g72
        assert 72*int(abs(rg).max())*int(abs(reps).max()) < 2**53
        maximum=-10**9
        for start in range(0,count,256):
            ip=rg[start:start+256].astype(float)@reps.T.astype(float)
            assert np.array_equal(ip,np.rint(ip))
            local=np.arange(len(ip));ip[local,start+local]=-10**9
            maximum=max(maximum,int(ip.max()))
        assert maximum<=6
        centre=np.stack([a[0],a[1],a[0]+a[1]+2*unbits(h)])[None,:,:]
        u=to_catalogue(centre,a,b)[0]
        assert not ((reps-u)%2).any()
        assert int(u@g72@u)%4==2
        stored['head_%d'%h]=reps.astype(np.int16)
        stored['centre_%d'%h]=u
        result[str(h)]={'nearest_subset_size':count,'representative_norm_squared':22,
                        'head_norm_squared':'11/2','largest_distinct_representative_inner_product':maximum}
    np.savez_compressed(HERE/'gamma72_coset_heads.npz',**stored)
    return result


def conference10():
    # GF(9)=GF(3)[t]/(t^2+1), encoded a+3b.
    def mul(x,y):
        a,b=x%3,x//3;c,d=y%3,y//3
        return (a*c-b*d)%3+3*((a*d+b*c)%3)
    def char(x):
        if x==0:return 0
        fourth=mul(mul(x,x),mul(x,x))
        assert fourth in (1,2)
        return 1 if fourth==1 else -1
    c=np.ones((10,10),dtype=np.int64);np.fill_diagonal(c,0)
    for x in range(9):
        for y in range(9):
            diff=(x%3-y%3)%3+3*((x//3-y//3)%3)
            c[x+1,y+1]=char(diff)
    assert np.array_equal(c,c.T) and np.array_equal(c@c,9*np.eye(10,dtype=np.int64))
    assert int(np.trace(c))==0
    return c


def tail_certificates():
    r=sp.sqrt
    v2=sp.Matrix([[1,0],[sp.Rational(1,2),r(3)/2],[-sp.Rational(1,2),r(3)/2]])
    v3=sp.Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/r(3)
    v4=sp.Matrix([[r(2)/3,sp.Rational(2,3)*sp.cos(2*sp.pi*i/3),
                   sp.Rational(2,3)*sp.sin(2*sp.pi*i/3),(-1)**i/r(3)] for i in range(6)])
    for v,dim,coherence in [(v2,2,sp.Rational(1,2)),(v3,3,sp.Rational(1,3)),(v4,4,sp.Rational(1,3))]:
        gram=sp.simplify(v*v.T)
        assert gram.rank()==dim and all(gram[i,i]==1 for i in range(len(gram.diagonal())))
        assert max(abs(gram[i,j]) for i in range(gram.rows) for j in range(i))<=coherence
    # Pure-tail roots D3 and D4 have pairwise unit products <=1/2.
    for v,dim in [(v3,3),(v4,4)]:
        maximum=sp.Integer(0)
        for i,j in itertools.combinations(range(dim),2):
            for si,sj in itertools.product((-1,1),repeat=2):
                z=sp.zeros(dim,1);z[i]=si/r(2);z[j]=sj/r(2)
                maximum=max(maximum,*[sp.simplify(abs(x)) for x in v*z])
        assert sp.simplify(maximum**2-sp.Rational(4,5))<0
    # Rotated hexagon: |cos(angle)|<=sqrt(3)/2 with the three tail lines.
    assert sp.Rational(3,4)<sp.Rational(4,5)
    c=conference10()
    # C^2=9I and tr C=0 certify I+C/3 is PSD of rank five.
    assert sp.Rational(5,1)+5*(1-sp.Rational(1,3))>=8
    return {'74':{'coset_count':3,'tail_coherence':'1/2','coset_distance_squared_at_least':'11/2','pure_tail_points':6},
            '75':{'coset_count':4,'tail_coherence':'1/3','coset_distance_squared_at_least':5,'pure_tail_points':12},
            '76':{'coset_count':6,'tail_coherence':'1/3','coset_distance_squared_at_least':5,'pure_tail_points':24},
            '77':{'coset_count':10,'tail_coherence':'1/3','coset_distance_squared_at_least':5,'pure_tail_points':0}}


def main():
    started=time.time();hashes={}
    structure=np.load(HERE/'gamma72_leech_structure.npz')
    f,a,b,t=[structure[k] for k in ('gram','a','b','t')]
    g72=np.load(SOURCE/'gamma72_gram.npy').astype(np.int64)
    assert np.array_equal(t@np.kron(np.eye(3,dtype=np.int64),f)@t.T,2*g72)
    assert np.array_equal(a@f@a.T,2*f)
    assert np.array_equal(a@b,2*np.eye(24,dtype=np.int64))
    assert sp.Matrix(f).det()==1 and not (np.diag(f)%2).any()
    check_intervals(f)
    full_g=np.load(SOURCE/'gamma72_gens.npy')[1].astype(np.int64)
    g=full_g[:24,:24]
    assert np.array_equal(full_g,np.kron(np.eye(3,dtype=np.int64),g))
    assert np.array_equal(full_g@g72@full_g.T,g72)
    assert np.array_equal(g@a,a@g)
    assert np.array_equal(np.linalg.matrix_power(g,12),-np.eye(24,dtype=np.int64))
    qs=[a[0],a[1],a[0]+a[1]]
    absence={}
    for power in range(7):
        mat=np.eye(24,dtype=np.int64) if not power else np.eye(24,dtype=np.int64)-np.linalg.matrix_power(g,power)
        qfamily=[q@mat for q in qs]
        assert all(np.any(q%2) and not (q@b%2).any() for q in qfamily)
        targets=[bits(unbits(h)@mat) for h in HEADS] if not power else [bits(unbits(6896114)@mat)]
        if power in (2,4,6):targets.append(bits(unbits(12625758)@mat))
        if power==4:
            triangle_h=bits(unbits(4915087)@mat);targets.append(triangle_h)
            actual=sp.Rational(sum(int(q@f@q) for q in qfamily),4)+int(qfamily[2]@f@unbits(triangle_h))
            assert actual%2==1
        records=load_lists(power,hashes,qfamily)
        absence[str(power)]=no_shorter_triples(records,a,sorted(set(targets)))
        print('verified absence of shorter vectors: difference power',power,flush=True)
    heads=check_heads(a,b,f,g72,hashes)
    tails=tail_certificates()
    specs={73:(12625758,1,0),74:(4915087,3,6),75:(12625758,4,12),
           76:(12625758,6,24),77:(6896114,10,0)}
    bounds={}
    for n,(h,m,poles) in specs.items():
        count=BASE+2*HEADS[h]*m+poles
        bounds[str(n)]={'bound':count,'repository_baseline':PREVIOUS[n],
                        'gain':count-PREVIOUS[n],'head_subset_size':HEADS[h],
                        'cosets':m,'pure_tail_points':poles,'h':h}
        assert count>PREVIOUS[n]
    for name in ['gamma72_coset_heads.npz','gamma72_leech_ldl.txt','gamma72_leech_structure.npz',
                 'gamma72_coset_certificate.py','gamma72_leech_enum.cpp','gamma72_coset_screen.cpp',
                 'gamma72_coset_count.cpp','gamma72_coset_extract.cpp']:
        hashes[name]=sha(HERE/name)
    report={'status':'Constructive lower bounds; completed interval enumeration plus independent direct checks',
            'base_lattice':'Gamma72, minimum squared norm 8','base_shell_size':BASE,
            'head_sets':heads,'shorter_vector_certificates':absence,'tail_codes':tails,
            'dimensions':bounds,'seconds':time.time()-started,'sha256':hashes,
            'scope':'Comparison is against the pinned repository benchmarks; no global priority claim.',
            'remaining':'The present certificates do not improve dimensions 78--94.'}
    (HERE/'gamma72_coset_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    for n,entry in bounds.items():print(n,format(entry['bound'],','),'gain',format(entry['gain'],','))


if __name__=='__main__':main()
