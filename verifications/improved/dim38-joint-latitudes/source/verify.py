"""Exact independent verifier for the accompanying 38D kissing construction.

Requires NumPy. Uses integer arithmetic exclusively for mathematical checks.
The reported elapsed time is the only floating-point output. No optimization
or numerical orthogonality tolerance is used. Run without Python's -O flag.
"""
from pathlib import Path
from itertools import combinations
from collections import Counter
import hashlib,json,time
import numpy as np
HERE=Path(__file__).resolve().parent
def pc(x):return bin(x).count('1')
def require(test,message):
    if not test:raise ValueError(message)
def multiply(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def transpose(A):return [list(c) for c in zip(*A)]
def leech():
    rows=(HERE/'source/golay_basis.txt').read_text().splitlines()
    require(len(rows)==12 and all(len(s)==24 and set(s)<={'0','1'} for s in rows),'Golay format')
    words=[0]
    for s in rows:
        b=sum((c=='1')<<i for i,c in enumerate(s));words += [w^b for w in words]
    require(len(set(words))==4096,'Golay rank')
    require(dict(Counter(map(pc,words)))=={0:1,8:759,12:2576,16:759,24:1},'Golay enumerator')
    vectors=[]
    for a,b in combinations(range(24),2):
        for sign in (-1,1):
            v=[0]*24;v[a]=4;v[b]=sign*4;vectors.append(tuple(v))
    for word in words:
        if pc(word)!=8:continue
        p=[i for i in range(24) if word>>i&1]
        for mask in range(128):
            if pc(mask)%2:continue
            v=[0]*24;v[p[0]]=2
            for j in range(7):v[p[j+1]]=2-4*((mask>>j)&1)
            vectors.append(tuple(v))
    for word in words:
        base=[1-2*((word>>i)&1) for i in range(24)]
        for a in range(24):
            v=base.copy();v[a]*=-3
            if v[0]>0:vectors.append(tuple(v))
    require(len(vectors)==len(set(vectors))==98280,'Leech line count/distinctness')
    vectors.sort();A=np.array(vectors,dtype=np.int16)
    require(np.all(np.sum(A.astype(np.int64)**2,axis=1)==32),'Leech norms')
    return A,set(words)
def membership(A,words):
    A=A.astype(np.int64);parity=A[:,0]%2
    require(np.all(A%2==parity[:,None]),'Leech common coordinate parity')
    residues=(((A-parity[:,None])//2)%2)@(1<<np.arange(24,dtype=np.int64))
    require(all(int(r) in words for r in residues),'Leech Golay residue')
    require(np.all(np.sum(A,axis=1)%8==4*parity),'Leech coordinate sum')
def surd_inequality(a,b,c):
    # a+b*sqrt(2) <= c*sqrt(6), with nonnegative a,b,c.
    remainder=6*c*c-a*a-2*b*b
    return remainder>=0 and remainder*remainder>=8*a*a*b*b
def main():
    if not __debug__:raise RuntimeError('Run without -O.')
    start=time.monotonic();cert=json.loads((HERE/'certificate.json').read_text())
    for name,expected in cert['source_sha256'].items():
        require(hashlib.sha256((HERE/'source'/name).read_bytes()).hexdigest()==expected,'Source hash: '+name)
    A,words=leech();membership(A,words)
    # The Golay enumerator and lattice congruences prove minimum norm 32;
    # see proof.txt for the even-coordinate and odd-coordinate cases.
    labels=np.loadtxt(HERE/'source/class_labels.txt',dtype=np.int64)
    W=np.loadtxt(HERE/'source/dim14_directions.txt',delimiter=',',dtype=np.int64)
    triples=np.loadtxt(HERE/'source/triples.txt',dtype=np.int64)
    require(labels.shape==(98280,) and np.all((labels>=0)&(labels<644)),'Class labels')
    require(W.shape==(1932,14) and np.all(np.sum(W*W,axis=1)==8),'Head norms')
    require(len(set(map(tuple,W.tolist())))==1932,'Head distinctness')
    G=W@W.T;np.fill_diagonal(G,-8);require(int(G.max())<=4,'Head pair angles')
    require(triples.shape==(644,3) and sorted(triples.ravel().tolist())==list(range(1932)),'Head partition')
    require(np.all(np.sum(W[triples],axis=1)==0),'Head triangle sums')
    max_class=0
    for c in range(644):
        B=A[labels==c];G=B@B.T;np.fill_diagonal(G,0)
        if len(B):max_class=max(max_class,int(np.abs(G).max()))
    require(max_class<=8,'Within-class body angles')
    print('Verified all 589680 original caps.',flush=True)

    E=np.loadtxt(HERE/'source/equator_lines_6855.txt',dtype=np.int64)
    require(E.shape==(6855,24) and np.all(np.sum(E*E,axis=1)==48),'Equator norms')
    membership(E,words);E=E.astype(np.int16)
    first=np.argmax(E!=0,axis=1);E[E[np.arange(len(E)),first]<0]*=-1
    require(len(set(map(tuple,E.tolist())))==6855,'Equator lines distinct')
    max_equator=0
    for begin in range(0,len(E),256):
        G=E[begin:begin+256]@E.T
        for i in range(len(G)):G[i,begin+i]=0
        max_equator=max(max_equator,int(np.abs(G).max()))
    require(max_equator<=24,'Equator angles')
    print('Verified all 13710 equatorial points; cap compatibility follows from lattice minimum.',flush=True)

    N=np.array(cert['rotation_numerator'],dtype=np.int64)
    require(cert['rotation_denominator']==1 and N.shape==(14,14),'Integer rotation format')
    require(np.all(np.abs(N)<=1) and np.array_equal(N@N.T,np.eye(14,dtype=np.int64)),'Exact orthogonality')
    keep=cert['axis_keep']
    require(len(keep)==len(set(keep)) and all(type(i) is int and 0<=i<len(W) for i in keep),'Axis selection')
    axis=W[keep]@N.T
    maximum_axis_dot=int(np.max(np.abs(W@axis.T)))
    require(maximum_axis_dot**2<=48,'Axis versus cap exact inequality')
    print('Verified',len(keep),'axis points by exact integer orthogonality and products.',flush=True)

    Z=np.array(cert['doubled_head_lines'],dtype=np.int64);body_ids=cert['body_line_indices']
    require(Z.ndim==2 and Z.shape[1]==14 and len(Z)==len(body_ids)>0,'New family dimensions')
    require(all(type(i) is int and 0<=i<98280 for i in body_ids),'Body indices')
    require(np.all(np.sum(Z*Z,axis=1)==24),'New head norms')
    U=A[body_ids].astype(np.int64)
    # Body lines may repeat. Check the joint body/head inequality directly.
    G=np.abs(U@U.T)+4*np.abs(Z@Z.T);np.fill_diagonal(G,0)
    maximum_new_pair=int(G.max());require(maximum_new_pair<=64,'Every pair of new families')
    require(int(np.abs(Z@W.T).max())<=8,'New heads versus old cap heads')
    maximum_new_axis=int(np.abs(Z@axis.T).max())
    require(maximum_new_axis<=8,'New heads versus retained axes')
    for i,z in zip(body_ids,Z):
        require(np.max(np.abs(W[triples[labels[i]]]@z))<=2,'New head versus its own body class')
    require(surd_inequality(16,16,16),'Different-body cap inequality')
    require(surd_inequality(32,4,16),'Same-body cap inequality')
    # Within each four-point family the distinct products are 64,-64,-128.
    # Against any equatorial point, Cauchy-Schwarz applies: body fraction 1/4.
    coordinate_axes=cert['additional_coordinate_axes']
    require(coordinate_axes==list(range(14)),'Additional coordinate axes')
    # Add (0,+/-e_j). Pair conditions reduce to coordinate-square bounds.
    require(np.max(W*W)<=6,'Coordinate axes versus caps')
    require(np.max(axis*axis)<=2,'Coordinate axes versus retained axes')
    require(np.max(Z*Z)<=8,'Coordinate axes versus new latitude')
    coordinate_count=2*len(coordinate_axes)
    total=6*len(A)+2*len(E)+len(keep)+coordinate_count+4*len(Z)
    require(cert['dimension']==38 and total==cert['claimed_size'] and total>cert['previous_size'],'Count')
    report={'all_checks_passed':True,'dimension':38,'lower_bound':total,
            'public_benchmark':cert['previous_size'],'gain_over_public':total-cert['previous_size'],
            'caps':6*len(A),'equator_points':2*len(E),'axis_points':len(keep)+coordinate_count,
            'retained_flat_axis_points':len(keep),'additional_coordinate_axis_points':coordinate_count,
            'new_points':4*len(Z),'new_families':len(Z),'distinct_new_body_lines':len(set(body_ids)),
            'max_within_class_inner_product':max_class,'max_equator_absolute_inner_product':max_equator,
            'max_axis_cap_head_product':maximum_axis_dot,'max_new_family_product':maximum_new_pair,
            'max_doubled_head_axis_product':maximum_new_axis,
            'arithmetic':'Integer arithmetic only; no numerical tolerances',
            'certificate_sha256':hashlib.sha256((HERE/'certificate.json').read_bytes()).hexdigest(),
            'elapsed_seconds':time.monotonic()-start}
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));print('VERIFIED: tau_38 >= '+str(total)+'.')
if __name__=='__main__':main()
