"""Independent expanded-coordinate proof; imports no symbolic/baseline verifier."""
import numpy as np,json,itertools,hashlib,time,argparse
from pathlib import Path
HERE=Path(__file__).resolve().parent

def verify(certificate=HERE/'certificate.json',write_report=False,coordinate_output=None):
    assert __debug__, 'Run without -O'
    p=Path(certificate).resolve().parent
    started=time.time();c=json.loads(Path(certificate).read_text());assert all(type(x)==int for key in ('stars','heavy_numerators_over_85') for row in c[key] for x in row);S=np.array(c['stars'],np.int64);H=np.array(c['heavy_numerators_over_85'],np.int64);coords=c['dense_coordinates'];assert coords==[j for j in range(34) if j not in [0,2]]
    assert c['dimension']==34 and c['points']==390992
    assert c['common_coordinate_scale']==425 and c['squared_norm']==5780000 and c['maximum_inner_product']==2890000 and c['heavy_multiplier_at_scale425']==4
    W=c['full_even_supports'];assert len(W)==len(set(W))==2013;assert all(type(w)==int and w>=0 and w<2**34 and w.bit_count()==8 for w in W)
    maxov=max((a&b).bit_count() for i,a in enumerate(W) for b in W[:i]);assert maxov<=4
    pat=np.array([s for s in itertools.product([-2,2],repeat=8) if sum(t<0 for t in s)%2==0],np.int64);assert pat.shape==(128,8)
    pg=pat@pat.T;np.fill_diagonal(pg,-999);assert pg.max()==16
    M=np.zeros((len(W)*128,34),np.int64)
    for i,m in enumerate(W):M[i*128:(i+1)*128,[j for j in range(34) if m>>j&1]]=pat
    words=[0]
    for g in c['dense_generator']:
     assert type(g)==int and 0<=g<2**32
     words += [w^g for w in words]
    assert len(words)==len(set(words))==131072
    mind=min(w.bit_count() for w in words[1:]);assert mind>=8
    shift=c['dense_coset_shift'];assert type(shift)==int and 0<=shift<2**32
    D=np.zeros((len(words),34),np.int64);D[:,coords]=1-2*(((np.array(words,np.uint32)^np.uint32(shift))[:,None]>>np.arange(32,dtype=np.uint32))&1).astype(np.int64)
    R=[]
    for a,b in itertools.combinations(range(34),2):
     for s,t in itertools.product([-4,4],repeat=2):
      v=np.zeros(34,np.int64);v[a]=s;v[b]=t;R.append(v)
    R=np.array(R);keep=(R@S.T).max(1)<=16;keep&=(R@H.T).max(1)<=1700;removed=R[~keep];R=R[keep]
    removed_lines=set()
    for row in removed:
     a,b=np.flatnonzero(row);removed_lines.add((int(a),int(b),int(row[a]*row[b]//16)))
    assert removed_lines==set(map(tuple,c['removed_root_lines'])) and len(removed_lines)*2==len(removed)
    assert np.all((H*H).sum(1)==361250)
    for layer in [D,M,R,S]:assert np.all((layer*layer).sum(1)==32)
    maxima={'dense_dense':32-2*mind,'middle_middle':max(4*maxov,int(pg.max()))}
    # Bound dense interactions coordinatewise. Every dense nonzero coordinate is +/-1.
    for name,L in [('middle',M),('roots',R),('stars',S)]:
     v=int(np.abs(L[:,coords]).sum(1).max());assert v<=16;maxima['dense_'+name+'_upper_bound']=v
    # Exact remaining ordinary mixed pairs and small same-layer Grams.
    for name,A,B in [('stars_middle',S,M),('stars_roots',S,R),('roots_middle',R,M)]:
     if name=='roots_middle':
      # Root has two entries +/-4, middle entries bounded by2: dot<=16.
      assert np.all((np.abs(R)==4).sum(1)==2) and np.abs(M).max()==2;v=16
     else:v=int((A@B.T).max())
     assert v<=16;maxima[name]=v
    for name,L in [('stars',S),('roots',R)]:
     g=L@L.T;np.fill_diagonal(g,-999);v=int(g.max());assert v<=16;maxima[name+'_'+name]=v
    # Heads are 4H/425; norm32. Mixed products<=16 iff H.L<=1700.
    for name,L in [('dense',D),('middle',M),('stars',S),('roots',R)]:
     v=int((H@L.T).max());assert v<=1700;maxima['head_numerator_'+name]=v
    hg=H@H.T;np.fill_diagonal(hg,-999999);assert hg.max()<=180625;maxima['head_numerator_head_numerator']=int(hg.max())
    # Common scale425 makes all coordinates integers, all norms5780000.
    ALL=np.concatenate([425*D,425*M,425*R,425*S,4*H]).astype(np.int32)
    assert len(ALL)==c['points']==390992
    assert c['counts']==dict(dense=len(D),middle=len(M),roots=len(R),stars=len(S),heavy=len(H))
    assert len(c['removed_root_lines'])==len(removed_lines)==28
    assert len(np.unique(ALL,axis=0))==len(ALL)
    assert np.all((ALL.astype(np.int64)**2).sum(1)==5780000)
    # Full rank follows from an exact rational row elimination on selected roots.
    from fractions import Fraction
    basis={}
    for row in R:
     v=list(map(Fraction,map(int,row)))
     for k,b in sorted(basis.items()):
      if v[k]:q=v[k];v=[a-q*z for a,z in zip(v,b)]
     k=next((i for i,a in enumerate(v) if a),None)
     if k is not None:
      q=v[k];basis[k]=[a/q for a in v]
     if len(basis)==34:break
    assert len(basis)==34
    out={'all_checks_pass':True,'points':len(ALL),'dimension':34,'exact_rank':len(basis),'counts':{'dense':len(D),'middle':len(M),'roots':len(R),'stars':len(S),'heads':len(H)},'removed_roots':len(removed),'dense_minimum_distance':mind,'middle_maximum_overlap':maxov,'maxima':maxima,'common_scale':425,'integer_squared_norm':5780000,'integer_pairwise_bound':2890000,'unique_rows':len(ALL),'input_sha256':{'certificate.json':hashlib.sha256(Path(certificate).read_bytes()).hexdigest()},'elapsed_seconds':time.time()-started,'inherited_assumptions':[],'scope_note':'Exact finite construction verified independently; current world-record/literature priority not established by this computation.'}
    if write_report:(HERE/'coordinate_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    if coordinate_output is not None:np.savez_compressed(coordinate_output,coordinates=ALL)
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--coordinates',help='Optional output .npz path for all integer coordinates')
    args=parser.parse_args()
    print(json.dumps(verify(write_report=True,coordinate_output=args.coordinates),indent=2))
