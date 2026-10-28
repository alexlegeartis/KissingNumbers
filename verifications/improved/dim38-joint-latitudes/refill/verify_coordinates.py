"""Independent expanded-coordinate proof; no imports from either other verifier."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import argparse
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'source/source'

def check(condition, message):
    if not condition:
        raise ValueError(message)

def reconstruct():
    rows=(SOURCE/'golay_basis.txt').read_text().splitlines()
    check(len(rows)==12 and all(len(s)==24 and set(s)<={'0','1'} for s in rows),'Golay basis format')
    generators=[sum(int(c)<<i for i,c in enumerate(s)) for s in rows]
    words=set()
    for mask in range(4096):
        word=0
        for i,g in enumerate(generators):
            if mask>>i&1:
                word^=g
        words.add(word)
    weight=lambda w: bin(w).count('1')
    check(len(words)==4096 and Counter(map(weight,words))=={0:1,8:759,12:2576,16:759,24:1},'Golay rank and enumerator')
    vectors=set()
    for a,b in combinations(range(24),2):
        for sign in (1,-1):
            x=[0]*24
            x[a],x[b]=4,4*sign
            vectors.add(tuple(x))
    for w in words:
        if weight(w)==8:
            support=[i for i in range(24) if w>>i&1]
            for signs in product((-1,1),repeat=7):
                if sum(s<0 for s in signs)%2:
                    continue
                x=[0]*24
                x[support[0]]=2
                for t,s in zip(support[1:],signs):
                    x[t]=2*s
                vectors.add(tuple(x))
        base=[1-2*((w>>i)&1) for i in range(24)]
        for t in range(24):
            x=base.copy()
            x[t]*=-3
            if x[0]>0:
                vectors.add(tuple(x))
    check(len(vectors)==98280,'Independent Leech line count')
    return np.array(sorted(vectors),dtype=np.int32),words

def member(X,words):
    X=X.astype(np.int64)
    parity=X[:,0]%2
    check(np.all(X%2==parity[:,None]),'Lattice parity')
    residues=(((X-parity[:,None])//2)%2)@(2**np.arange(24,dtype=np.int64))
    check(all(int(r) in words for r in residues),'Lattice Golay residue')
    check(np.all(X.sum(axis=1)%8==4*parity),'Lattice coordinate sum')
    # These congruences and the enumerator give minimum norm 32:
    # nonzero even residue uses >=8 coordinates of magnitude >=2;
    # zero even residue uses multiples of 4, with a single +/-4 forbidden;
    # an odd vector below norm32 would have all magnitudes1 and wrong sum mod8.

def max_products(X,absolute):
    maximum=-10**9
    for start in range(0,len(X),128):
        G=X[start:start+128]@X.T
        if absolute:
            G=np.abs(G)
        for j in range(len(G)):
            G[j,start+j]=-10**9
        maximum=max(maximum,int(G.max()))
    return maximum

def main():
    check(__debug__,'Run without -O')
    parser=argparse.ArgumentParser()
    parser.add_argument('--certificate',type=Path,default=HERE/'certificate.json')
    args=parser.parse_args()
    c=json.loads(args.certificate.read_text())
    check(c.get('format')=='fano_grid_cube_latitudes_v1','Certificate format')
    baseline=json.loads((HERE.parent/'source/certificate.json').read_text())
    for field in ('dimension','source_sha256','rotation_numerator','rotation_denominator','axis_keep','additional_coordinate_axes'):
        check(c[field]==baseline[field],'Retained field: '+field)
    for name,digest in c['source_sha256'].items():
        check(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==digest,'Source hash: '+name)
    A,words=reconstruct()
    member(A,words)
    check(np.all(np.sum(A*A,axis=1)==32),'Body norms')
    labels=np.loadtxt(SOURCE/'class_labels.txt',dtype=np.int32)
    W=np.loadtxt(SOURCE/'dim14_directions.txt',delimiter=',',dtype=np.int32)
    triples=np.loadtxt(SOURCE/'triples.txt',dtype=np.int32)
    check(labels.shape==(98280,) and np.all((labels>=0)&(labels<644)),'Cap labels')
    check(W.shape==(1932,14) and np.max(np.abs(W))<=2 and np.all((W*W).sum(axis=1)==8),'Cap head norms')
    check(len({tuple(x) for x in W})==1932 and max_products(W,False)<=4,'Cap head distinctness and products')
    check(triples.shape==(644,3) and sorted(triples.ravel())==list(range(1932)),'Cap triangles')
    check(np.all(W[triples].sum(axis=1)==0),'Zero-sum cap triangles')
    for label in range(644):
        B=A[labels==label]
        check(len(B)<2 or max_products(B,True)<=8,'Within-class body products')
    E=np.loadtxt(SOURCE/'equator_lines_6855.txt',dtype=np.int32)
    check(E.shape==(6855,24) and np.max(np.abs(E))<=6 and np.all((E*E).sum(axis=1)==48),'Equator norms')
    member(E,words)
    canonical=[]
    for e in E:
        sign=1 if next(int(x) for x in e if x)>0 else -1
        canonical.append(tuple(sign*e))
    check(len(set(canonical))==6855 and max_products(E,True)<=24,'Equator distinctness and products')
    N=np.array(c['rotation_numerator'],dtype=np.int32)
    check(c['rotation_denominator']==1 and N.shape==(14,14) and np.max(np.abs(N))<=1,'Axis rotation format')
    check(np.array_equal(N@N.T,np.eye(14,dtype=np.int32)),'Axis orthogonality')
    keep=c['axis_keep']
    check(len(keep)==len(set(keep)) and all(type(i) is int and 0<=i<1932 for i in keep),'Axis indices')
    axis=W[keep]@N.T
    check(int(np.abs(W@axis.T).max())**2<=48,'Axis versus cap')
    check(c['additional_coordinate_axes']==list(range(14)),'Coordinate axes')
    check(np.max(W*W)<=6 and np.max(axis*axis)<=2,'Coordinate axes versus old points')

    ids=c['body_line_indices']
    heads=c['doubled_head_lines']
    check(len(ids)==len(heads)>0 and all(type(i) is int and 0<=i<98280 for i in ids),'Grid body format')
    removed=c['removed_cap_body_lines']
    check(len(removed)==len(set(removed)) and all(type(i) is int and 0<=i<98280 for i in removed),'Deleted cap owners')
    rows=[]
    owners=set(ids)
    for body,z in zip(ids,heads):
        check(len(z)==14 and all(type(v) is int and abs(v)<=2 for v in z) and sum(v*v for v in z)==24,'Grid head norm')
        for sb,sh in product((-1,1),repeat=2):
            rows.append([sb*int(v) for v in A[body]]+[2*sh*v for v in z])
    extras=c['additional_latitude_points']
    for e in extras:
        body,sign,z=e['body_line_index'],e['body_sign'],e['doubled_head']
        check(type(body) is int and 0<=body<98280 and type(sign) is int and sign in (-1,1),'Individual body format')
        check(len(z)==14 and all(type(v) is int and abs(v)<=2 for v in z) and sum(v*v for v in z)==24,'Individual head norm')
        owners.add(body)
        rows.append([sign*int(v) for v in A[body]]+[2*v for v in z])
    check(set(removed)==owners,'Every latitude body must lose all six same-body caps')
    P=np.array(rows,dtype=np.int32)
    check(P.shape==(4*len(ids)+len(extras),38) and np.max(np.abs(P))<=4,'Expanded point format')
    check(np.all((P*P).sum(axis=1)==128) and np.all((P[:,:24]*P[:,:24]).sum(axis=1)==32),'Expanded norms')
    check(len({tuple(row) for row in P})==len(P),'Expanded point distinctness')
    # Every entry is <=4, so each product accumulator is <=38*16=608,
    # well inside signed int32; all matrix arithmetic below is exact.
    maximum=max_products(P,False)
    check(maximum<=64,'All distinct expanded latitude pairs')
    tail=P[:,24:]
    check(int(np.abs(tail@W.T).max())<=16,'Latitude versus retained cap heads')
    check(int(np.abs(tail@axis.T).max())<=16,'Latitude versus flat axes')
    check(np.max(tail*tail)<=32,'Latitude versus coordinate axes')
    # For any retained cap, its line is different from the latitude body.
    # Leech minimum gives |u.u_cap|<=16. Thus the numerator is at most
    # 16+16sqrt(2) <= 16sqrt(6), with denominator32sqrt(6).
    # Prove this independently by sqrt(2)<=10/7 and sqrt(6)>=17/7:
    check(100>=2*49 and 17*17<=6*49 and 7*16+10*16<=17*16,'Exact different-body cap bound')
    # Equatorial products are <=1/2 by Cauchy-Schwarz on a body of
    # length1/2. The cap/equator bound follows from minimum32 and norm48.
    total=6*(len(A)-len(removed))+2*len(E)+len(axis)+28+len(P)
    check(c['dimension']==38 and total==c['claimed_size'] and total>c['previous_local_size'],'Count')
    report=dict(all_checks_passed=True,dimension=38,lower_bound=total,
        retained_caps=6*(len(A)-len(removed)),removed_caps=6*len(removed),
        latitude_points=len(P),latitude_pairs=len(P)*(len(P)-1)//2,
        max_latitude_inner_product_numerator=maximum,
        certificate_sha256=hashlib.sha256(args.certificate.read_bytes()).hexdigest(),
        arithmetic='Exact int32 products, int64 lattice congruences; no tolerances')
    (HERE/'coordinate_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    print('VERIFIED independently: tau_38 >= '+str(total))

if __name__=='__main__':
    main()
