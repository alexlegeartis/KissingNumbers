"""Independent NumPy integer-coordinate verifier of the complete frozen code.
No search, graph or scalar verifier is imported. Only entries, integer
products and exact shape bounds are used. Full all-pairs Gram is unnecessary.
"""
from pathlib import Path
import hashlib,itertools,json,sys,time
import numpy as np

def verify(path):
 start=time.monotonic();c=json.loads(path.read_text());assert c['dimension']==33 and c['format']=='signed_octad_star_points_explicit_v1'
 assert c['dense_coordinates']==list(range(32)) and c['dense_coset_shift']==0
 words=np.array([0],dtype=np.uint32)
 for g in c['dense_generator']:words=np.concatenate((words,words^np.uint32(g)))
 assert len(words)==len(np.unique(words))==131072
 bits=((words[:,None]>>np.arange(32,dtype=np.uint32))&1).astype(np.int64)
 weights=bits.sum(axis=1);assert np.min(weights[weights!=0])>=8
 dense=np.zeros((len(words),33),dtype=np.int64);dense[:,:32]=1-2*bits
 patterns=np.array([s for s in itertools.product((-2,2),repeat=8) if sum(x<0 for x in s)%2==0],dtype=np.int64)
 gram=patterns@patterns.T;np.fill_diagonal(gram,-32);assert gram.max()==16
 full=[];incidence=[]
 for m in c['full_even_supports']:
  js=np.array([j for j in range(33) if m>>j&1]);assert len(js)==8 and 0<=m<1<<33
  layer=np.zeros((128,33),dtype=np.int64);layer[:,js]=patterns;full.append(layer)
  incidence.append(np.abs(layer[0])//2)
 F=np.array(incidence,dtype=np.int64);overlaps=F@F.T;np.fill_diagonal(overlaps,0);assert overlaps.max()<=4
 full=np.concatenate(full)
 E=[]
 for m,s in c['exceptional_signed_points']:
  assert 0<m<1<<33 and 0<=s and s&m==s
  E.append([(-2 if s>>j&1 else 2) if m>>j&1 else 0 for j in range(33)])
 E=np.array(E,dtype=np.int64);assert np.all(np.sum(E*E,axis=1)==32)
 assert np.max((np.abs(E)//2)@F.T)<=4
 eg=E@E.T;np.fill_diagonal(eg,-32);assert eg.max()<=16
 removed=set(map(tuple,c['removed_root_signed_points']));assert len(removed)==len(c['removed_root_signed_points'])
 assert all(len(v)==33 and sorted(abs(x) for x in v if x)==[4,4] for v in removed)
 roots=[];seenremoved=set()
 for a,b in itertools.combinations(range(33),2):
  for sa,sb in itertools.product((-4,4),repeat=2):
   v=[0]*33;v[a],v[b]=sa,sb
   if tuple(v) not in removed:roots.append(v)
   else:seenremoved.add(tuple(v))
 assert seenremoved==removed
 roots=np.array(roots,dtype=np.int64);rg=roots@roots.T;np.fill_diagonal(rg,-32);assert rg.max()<=16
 H=c['star_signed_points']
 assert all(len(v)==33 and sorted(abs(x) for x in v if x)==[2,2,2,2,4] for v in H)
 H=np.array(H,dtype=np.int64)
 middle=np.concatenate((full,E));maxima={};products=eg.size+rg.size
 for label,X in [('dense',dense),('full',full),('exceptional',E),('roots',roots)]:
  dot=H@X.T;assert dot.max()<=16;maxima[label]=int(dot.max());products+=dot.size
 hg=H@H.T;np.fill_diagonal(hg,-32);assert hg.max()<=16;maxima['stars']=int(hg.max());products+=hg.size
 arrays=[dense,full,E,roots,H]
 for X in arrays:assert np.all(np.sum(X*X,axis=1)==32)
 # Octads have eight entries +/-2; roots have two entries +/-4.
 # These shapes give middle/root <=16, middle/dense <=16, root/dense <=8.
 assert np.all(np.sum(middle!=0,axis=1)==8) and np.max(np.abs(middle))==2
 assert np.all(np.sum(roots!=0,axis=1)==2) and np.max(np.abs(roots))==4
 assert np.all(np.sum(dense!=0,axis=1)==32) and np.max(np.abs(dense))==1
 # Shape categories are disjoint; check uniqueness within each separately.
 for X in arrays:assert len(np.unique(X,axis=0))==len(X)
 assert not set(c['full_even_supports'])&{m for m,s in c['exceptional_signed_points']}
 total=sum(len(X) for X in arrays);assert total==c['points']
 # Exact mod-101 rank lower bound for real span, using retained roots.
 rows={}
 for v in roots:
  row=[int(x)%101 for x in v]
  for j in range(33):
   if row[j]:
    if j in rows:
     q=row[j];row=[(x-q*y)%101 for x,y in zip(row,rows[j])]
    else:
     inv=pow(row[j],-1,101);rows[j]=[(x*inv)%101 for x in row];break
  if len(rows)==33:break
 assert len(rows)==33
 return {'all_checks_pass':True,'points':total,'counts':dict(zip(['dense','full_octads','exceptional_points','roots','stars'],map(len,arrays))),'coordinate_dtype':'int64','explicit_integer_products':products,'star_product_maxima':maxima,'exceptional_pair_maximum':int(eg.max()),'root_pair_maximum':int(rg.max()),'full_support_overlap_maximum':int(overlaps.max()),'rank_mod_101':len(rows),'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
if __name__=='__main__':
 if not __debug__:raise SystemExit('Run without -O.')
 path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('certificate33_370408_explicit.json')
 print(json.dumps(verify(path),indent=2))
