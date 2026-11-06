"""Exact symbolic verification of the rationally rotated 390992-point code."""
import hashlib,itertools,json,time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent

def verify(path=HERE/'certificate.json',write_report=False):
    assert __debug__, 'Run without -O'
    started=time.monotonic();path=Path(path);c=json.loads(path.read_text())
    assert c['dimension']==34 and c['points']==390992
    assert c['common_coordinate_scale']==425 and c['squared_norm']==5780000
    assert c['maximum_inner_product']==2890000 and c['heavy_multiplier_at_scale425']==4
    W=c['full_even_supports'];assert len(W)==len(set(W))==2013
    assert all(type(w)==int and 0<=w<2**34 and w.bit_count()==8 for w in W)
    maximum_overlap=0
    for i,w in enumerate(W):
        for v in W[:i]:
            q=(w&v).bit_count();assert q<=4;maximum_overlap=max(maximum_overlap,q)
    coords=c['dense_coordinates'];assert coords==[j for j in range(34) if j not in [0,2]]
    words=[0]
    for g in c['dense_generator']:
        assert type(g)==int and 0<=g<2**32
        words += [w^g for w in words]
    assert len(words)==len(set(words))==131072
    mind=min(w.bit_count() for w in words[1:]);assert mind>=8
    shift=c['dense_coset_shift'];assert type(shift)==int and 0<=shift<2**32
    shifted=[(w^shift).bit_count() for w in words];assert min(shifted)>=6 and max(shifted)<=26
    assert all(type(x)==int for key in ('stars','heavy_numerators_over_85') for row in c[key] for x in row)
    S=np.array(c['stars'],np.int64);H=np.array(c['heavy_numerators_over_85'],np.int64)
    assert S.shape==(64,34) and H.shape==(4,34)
    assert len(np.unique(S,axis=0))==64 and len(np.unique(H,axis=0))==4
    assert np.all((S*S).sum(1)==32) and np.all((H*H).sum(1)==361250)
    assert np.max(S@S.T-32*np.eye(64,dtype=np.int64))<=16
    assert np.max(H@H.T-361250*np.eye(4,dtype=np.int64))<=180625
    removed=set(map(tuple,c['removed_root_lines']))
    assert len(removed)==len(c['removed_root_lines'])==28
    assert all(0<=a<b<34 and s in (-1,1) for a,b,s in removed)
    R=[]
    for a,b in itertools.combinations(range(34),2):
        for s,t in itertools.product((-4,4),repeat=2):
            if (a,b,s*t//16) in removed:continue
            v=np.zeros(34,np.int64);v[a]=s;v[b]=t;R.append(v)
    R=np.array(R);assert R.shape==(2188,34)
    inc=np.array([[(w>>j)&1 for j in range(34)] for w in W],np.int64)
    star_middle_bound=int((2*np.abs(S)@inc.T).max());assert star_middle_bound<=16
    assert np.max(S@R.T)<=16 and np.max(H@S.T)<=1700 and np.max(H@R.T)<=1700
    head_middle=[]
    for h in H:
        maxima=[]
        for w in W:
            z=h[[j for j in range(34) if w>>j&1]]
            # Optimize over all 128 even sign patterns exactly.
            value=2*int(np.abs(z).sum())
            if np.count_nonzero(z<0)%2:value-=4*int(np.abs(z).min())
            assert value<=1700;maxima.append(value)
        head_middle.append(max(maxima))
        assert len(set(map(int,h[coords])))==1 and abs(int(h[coords[0]]))==85
        assert max(abs((32-2*k)*int(h[coords[0]])) for k in shifted)<=1700
    # Dense/middle, dense/root and dense/star follow coordinatewise bounds.
    assert np.abs(S[:,coords]).sum(1).max()<=16 and np.abs(R[:,coords]).sum(1).max()<=16
    # Roots have two entries +/-4, hence root/root and root/middle products<=16.
    # Each even middle fiber has Hamming distance>=2; other supports overlap<=4.
    counts=dict(dense=len(words),middle=128*len(W),stars=len(S),heavy=len(H),roots=len(R))
    assert counts==c['counts'] and sum(counts.values())==c['points']==390992
    # Disjointness: dense has32 nonzeros of magnitude1; middle has eight of2;
    # roots have two of4; stars are unique, norm32, and below the dense/root/middle
    # pair bounds. Heads have common-scale core entries340, unlike all other layers.
    assert np.all(np.count_nonzero(S,axis=1)==5)
    assert np.all(np.sort(np.abs(S),axis=1)[:,-5:]==np.array([2,2,2,2,4]))
    assert np.all(np.abs(4*H[:,coords])==340)
    out=dict(all_checks_pass=True,dimension=34,points=390992,counts=counts,
             support_pairs=len(W)*(len(W)-1)//2,maximum_support_overlap=maximum_overlap,
             dense_minimum_distance=mind,dense_shifted_weight_range=[min(shifted),max(shifted)],
             head_middle_numerator_maxima=head_middle,star_middle_upper_bound=star_middle_bound,
             integer_squared_norm=5780000,integer_pair_bound=2890000,
             certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             elapsed_seconds=time.monotonic()-started)
    if write_report:(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':print(json.dumps(verify(write_report=True),indent=2))
