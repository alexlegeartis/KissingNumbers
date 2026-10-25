"""Portable exact verifier for midpoint replacement kissing configurations.

Only Python and NumPy are needed. Rebuilds the Leech shell, loads explicit
coordinates, and checks all new/new pairs with integer arithmetic. Retained
bulk pairs use the standard Leech construction; bulk/new pairs use exact
midpoint convexity, with endpoints verified to be deleted Leech vectors.
The pure hexagon (26D) or cuboctahedron (27D) is justified by Cauchy-Schwarz.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import numpy as np


def leech_shell():
    generator=np.array([1,0,1,0,1,1,1,0,0,0,1,1],dtype=np.uint8)
    matrix=np.zeros((12,23),dtype=np.uint8)
    for i in range(12):matrix[i,i:i+12]=generator
    bits=((np.arange(4096)[:,None]>>np.arange(12))&1).astype(np.uint8)
    golay=np.zeros((4096,24),dtype=np.uint8)
    golay[:,:23]=(bits@matrix)%2;golay[:,23]=golay[:,:23].sum(1)%2
    assert dict(zip(*np.unique(golay.sum(1),return_counts=True)))=={0:1,8:759,12:2576,16:759,24:1}
    rows=[]
    for i,j in itertools.combinations(range(24),2):
        for a,b in itertools.product([4,-4],repeat=2):
            row=np.zeros(24,dtype=np.int8);row[i]=a;row[j]=b;rows.append(row)
    signs=np.array([[(-2 if (k>>i)&1 else 2) for i in range(8)]
                    for k in range(256) if bin(k).count('1')%2==0],dtype=np.int8)
    for word in golay[golay.sum(1)==8]:
        block=np.zeros((128,24),dtype=np.int8);block[:,np.flatnonzero(word)]=signs
        rows.extend(block)
    sign_words=np.where(golay, -1,1).astype(np.int8)
    for j in range(24):
        block=sign_words.copy();block[:,j]*=-3;rows.extend(block)
    m=np.array(rows,dtype=np.int8)
    assert m.shape==(196560,24)
    assert np.all(np.sum(m.astype(np.int32)**2,axis=1)==32)
    assert len({row.tobytes() for row in m})==len(m)
    return m


def verify(path):
    start=time.time();data=np.load(path,allow_pickle=False)
    new=data['new_points'].astype(np.int16)
    removed=data['removed_leech'].astype(np.int16)
    endpoints=data['midpoint_owner_indices'].astype(np.int64)
    d=new.shape[1];assert d in [26,27]
    m=leech_shell();shell={row.tobytes() for row in m}
    assert len(np.unique(removed,axis=0))==len(removed)
    assert all(row.astype(np.int8).tobytes() in shell for row in removed)
    assert endpoints.shape==(len(new),2) and np.all(endpoints>=0) and np.all(endpoints<len(removed))
    assert np.all(endpoints[:,0]!=endpoints[:,1])
    assert len(np.unique(new,axis=0))==len(new)
    assert np.all(np.sum(new.astype(np.int32)**2,axis=1)==32)
    assert np.all(np.sum(new[:,24:].astype(np.int32)**2,axis=1)==8)
    assert np.all(2*new[:,:24]==removed[endpoints].sum(axis=1))
    # Consequently every new foot clears every retained bulk point exactly:
    # both endpoint products are <=16, so their average is <=16.
    worst=-100;contacts=0
    for at in range(0,len(new),256):
        block=new[at:at+256]
        dot=block@new.T
        # Each norm is sqrt(32). The sum of absolute coordinate products is
        # <=32 by Cauchy-Schwarz, so no int16 partial sum can overflow.
        dot[np.arange(len(block)),at+np.arange(len(block))]=-100
        assert np.max(dot)<=16,(at,int(np.max(dot)))
        worst=max(worst,int(dot.max()));contacts+=int(np.sum(dot==16))
    # Every tail has norm sqrt(8); a pure pole has norm sqrt(32), so their
    # product is <=16. Six regular-hexagon poles, or twelve cuboctahedron
    # poles, themselves form a kissing configuration in the tail space.
    pure=6 if d==26 else 12
    total=196560-len(removed)+len(new)+pure
    result=dict(dimension=d,points=total,retained_leech=196560-len(removed),
                removed_leech=len(removed),new_points=len(new),pure_points=pure,
                maximum_new_new_inner_product=worst,new_new_contacts=contacts//2,
                certificate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                arithmetic='exact integer coordinates, squared norm 32, threshold 16',
                seconds=time.time()-start)
    print(json.dumps(result,indent=2))
    print('PASS: EXACT LOWER BOUND tau_%d >= %d'%(d,total))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('certificate',type=Path)
    args=parser.parse_args();verify(args.certificate)
