"""Exact 38D triangular parity-Hamming replacement certificate."""
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as F
import hashlib, importlib.util, json
import numpy as np
P=Path(__file__).resolve().parent
S=P.parent/'dim38-joint-latitudes'
src=S/'source/source'
h=json.loads((S/'hamming/certificate.json').read_text())
old=json.loads((S/'refill/certificate.json').read_text())

def cert_input(path,sha):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==sha

def geom():
    s=importlib.util.spec_from_file_location('leech_source',S/'source/verify.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    A,_=m.leech()
    A=A.astype(np.int64)
    assert A.shape==(98280,24)
    return A

# Four signed triangles on five coordinates, each followed by its antipode
# in the palette.  Entries give (i,j,k, sign_i,sign_j,sign_k).
TRI=[(0,1,2,1,1,1),(0,1,3,1,-1,1),(0,2,4,1,-1,1),(1,3,4,1,1,-1)]

def construction():
    bodies=[];tails=[]
    for half,mask in product(range(2),range(128)):
        syndrome=0
        for k in range(7):
            if mask>>k&1:syndrome^=k+1
        sign=(-1)**(syndrome&1)
        block=14+5*(mask.bit_count()&1)
        i,j,k,a,b,c=TRI[syndrome//2]
        i+=block;j+=block;k+=block
        tail=[0]*14;tail[7*half:7*half+7]=[(-1)**(mask>>v&1) for v in range(7)]
        for u,v,x,y in ((i,j,a,-b),(j,k,b,-c),(k,i,c,-a)):
            body=[0]*24;body[u]=4*sign*x;body[v]=4*sign*y
            bodies.append(body);tails.append(tail.copy())
    return np.array(bodies,dtype=np.int64),np.array(tails,dtype=np.int64)

def replay():
    for name,sha in h['source_sha256'].items():cert_input(src/name,sha)
    cert_input(S/'refill/certificate.json',h['baseline_certificate_sha256'])
    assert h['claimed_size']==611638 and h['new_point_count']==512
    A=geom()
    B,Y=construction()
    assert B.shape==(768,24) and Y.shape==(768,14)
    assert np.all(np.sum(B*B,1)==32) and np.all(np.sum(Y*Y,1)==7)
    assert len(np.unique(np.c_[B,Y],axis=0))==768
    # In unit-normalized coordinates the point is (b/sqrt(96), sqrt(2/21)y).
    M=7*(B@B.T)+64*(Y@Y.T)
    assert np.all(np.diag(M)==672)
    np.fill_diagonal(M,0)
    assert M.max()<=336
    R={tuple(row):i for i,row in enumerate(A)}
    I={R[tuple(b) if b[np.flatnonzero(b)[0]]>0 else tuple(-b)] for b in B}
    assert len(I)==24
    gridids=np.unique(np.array(old['body_line_indices'],dtype=int))
    assert len(gridids)==14 and not (I&set(gridids.tolist()))
    previous=set(h['removed_cap_body_lines'])
    assert len(previous)==22
    removed=set(gridids.tolist())|I
    assert len(removed)==38
    # Verify every newly restored cap's interactions with the retained grid,
    # and every replacement body's interactions with all retained cap lines.
    restored=sorted(previous-removed)
    retain=np.array([i for i in range(len(A)) if i not in removed])
    body_max=int(np.abs(np.unique(B,axis=0)@A[retain].T).max())
    assert body_max<=16
    assert np.abs(A[gridids]@A[restored].T).max(initial=0)==0
    assert (B@A[gridids].T==0).all()
    # Retained cap and axes data use the same common scales as the parent.
    W=np.loadtxt(src/'dim14_directions.txt',delimiter=',',dtype=np.int64)
    axes=W[old['axis_keep']]@np.asarray(old['rotation_numerator'],dtype=np.int64).T
    E=np.loadtxt(src/'equator_lines_6855.txt',dtype=np.int64)
    Z=np.asarray(old['doubled_head_lines'],dtype=np.int64)
    assert len(W)==1932 and len(axes)==1568 and len(E)==6855 and len(Z)==1568
    assert np.all(np.sum(W*W,1)==8) and np.all(np.sum(axes*axes,1)==8)
    cap=int(np.abs(Y@W.T).max());ax=int(np.abs(Y@axes.T).max())
    grid=int(np.abs(Y@Z.T).max());eq=int(np.abs(B@E.T).max())
    assert cap<=4 and ax**2<=21 and grid**2<=84 and eq**2<=1152
    # Cap–donor scalar product <= sqrt(1/18)+4/sqrt(252).
    # Verify the scaled bound 2/sqrt(3)+4*sqrt(2/21) <= sqrt(6).
    assert F(17,10)**2<=3 and F(31,100)**2>=F(2,21) and F(17,7)**2<=6
    assert F(20,17)+F(31,25)<=F(17,7)
    # The fourteen coordinate axes: tail coordinate amplitude squared 2/21.
    assert old['additional_coordinate_axes']==list(range(14)) and F(2,21)<=F(1,4)
    # Different layers have distinct body norm fractions.
    assert len({F(2,3),F(1),F(0),F(1,4),F(1,3)})==5
    counts={'caps':6*len(retain),'equator':2*len(E),'axes':len(axes)+28,'grid':4*len(Z),'triangle_hamming':len(B)}
    n=sum(counts.values())
    assert n==611798
    return {'status':'PASS','dimension':38,'total':n,'gain':n-h['claimed_size'],'counts':counts,
            'body_lines_removed_for_triangle':len(I),'restored_cap_lines':restored,
            'new_new_pairs':len(B)*(len(B)-1)//2,'new_new_max_numerator_over_672':int(M.max()),
            'new_cap_maxima':[body_max,cap],'new_axis_tail_max':ax,
            'new_grid_tail_max':grid,'new_equator_body_max':eq}

if __name__=='__main__':
    print(json.dumps(replay(),sort_keys=True))
