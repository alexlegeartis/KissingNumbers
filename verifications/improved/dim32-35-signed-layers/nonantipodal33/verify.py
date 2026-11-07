"""Self-contained exact verifier for 33D explicit signed octads plus root stars.
Python standard library only. No search, graph, matching, source certificate,
or previously recorded verification is trusted. No antipodal completion.
"""
from pathlib import Path
import collections,hashlib,itertools,json,sys,time

def verify(path):
 start=time.monotonic();raw=path.read_bytes();c=json.loads(raw)
 assert c['format']=='signed_octad_star_points_explicit_v1' and c['dimension']==33
 assert c['dense_coordinates']==list(range(32))
 words={0}
 for g in c['dense_generator']:
  assert type(g) is int and 0<g<1<<32
  old=len(words);words|={w^g for w in words};assert len(words)==2*old
 assert len(words)==131072
 distance=min(w.bit_count() for w in words if w);assert distance>=8
 assert 0<=c['dense_coset_shift']<1<<32
 full=c['full_even_supports'];fullset=set(full)
 assert len(fullset)==len(full)
 assert all(type(m) is int and 0<m<1<<33 and m.bit_count()==8 for m in full)
 exceptional=[tuple(x) for x in c['exceptional_signed_points']];eset=set(exceptional)
 assert len(eset)==len(exceptional)
 for m,s in exceptional:assert type(m) is int and type(s) is int and 0<m<1<<33 and m.bit_count()==8 and s>=0 and s&m==s and m not in fullset
 exsupports={m for m,s in exceptional}
 ff=max(((a&b).bit_count() for a,b in itertools.combinations(full,2)),default=0);assert ff<=4
 fe=max(((a&b).bit_count() for a in full for b in exsupports),default=0);assert fe<=4
 ee=-32
 for (a,s),(b,t) in itertools.combinations(exceptional,2):
  overlap=a&b;dot=4*(overlap.bit_count()-2*((s^t)&overlap).bit_count());ee=max(ee,dot);assert dot<=16
 roots=[];removed={tuple(r) for r in c['removed_root_signed_points']}
 assert len(removed)==len(c['removed_root_signed_points'])
 assert all(len(r)==33 and sum(x!=0 for x in r)==2 and all(x in (-4,0,4) for x in r) and sum(x*x for x in r)==32 for r in removed)
 seenremoved=set()
 for a,b in itertools.combinations(range(33),2):
  for sa,sb in itertools.product((-4,4),repeat=2):
   v=[0]*33;v[a],v[b]=sa,sb
   if tuple(v) not in removed:roots.append((a,b,sa,sb))
   else:seenremoved.add(tuple(v))
 assert seenremoved==removed
 stars=[tuple(v) for v in c['star_signed_points']]
 assert all(len(v)==33 and sorted(abs(x) for x in v if x)!=[ ] and sorted(abs(x) for x in v if x)==[2,2,2,2,4] for v in stars)
 assert len(stars)==len(set(stars)) and all(sum(x*x for x in v)==32 for v in stars)
 star_max=dict(full=-32,exceptional=-32,roots=-32,stars=-32,dense_l1_bound=0)
 for v in stars:
  star_max['dense_l1_bound']=max(star_max['dense_l1_bound'],sum(abs(v[j]) for j in c['dense_coordinates']))
  for m in full:
   entries=[v[j] for j in range(33) if m>>j&1]
   score=2*sum(abs(x) for x in entries)
   if sum(x<0 for x in entries)%2:score-=4*min(abs(x) for x in entries)
   star_max['full']=max(star_max['full'],score);assert score<=16
  for m,s in exceptional:
   dot=2*sum(v[j]*(-1 if s>>j&1 else 1) for j in range(33) if m>>j&1)
   star_max['exceptional']=max(star_max['exceptional'],dot);assert dot<=16
  for a,b,sa,sb in roots:
   dot=sa*v[a]+sb*v[b];star_max['roots']=max(star_max['roots'],dot);assert dot<=16
 for v,w in itertools.combinations(stars,2):
  dot=sum(x*y for x,y in zip(v,w));star_max['stars']=max(star_max['stars'],dot);assert dot<=16
 assert star_max['dense_l1_bound']<=16
 counts={'dense':len(words),'full_even_octad_supports':len(full),'full_even_octad_points':128*len(full),'individual_exceptional_points':len(exceptional),'roots':len(roots),'stars':len(stars)}
 total=len(words)+128*len(full)+len(exceptional)+len(roots)+len(stars);assert total==c['points']
 # Remaining category bounds follow from shapes, norms and distinctness:
 # roots/roots <=16; root/octad <=16; dense/root <=8; dense/octad <=16.
 # Full octads on one support have even sign parity, so distinct points
 # differ in at least two positions and have product <=16.
 maxima={'dense_dense':32-2*distance,'full_same_support':16,'full_full':4*ff,'full_exceptional':4*fe,'exceptional_exceptional':ee,'roots_roots':16,'roots_octad':16,'dense_roots':8,'dense_octad':16,'star_categories':star_max}
 return {'all_checks_pass':True,'dimension':33,'points':total,'counts':counts,'unpaired_points':sum((m,s^m) not in eset for m,s in exceptional)+sum(tuple(-x for x in v) not in set(stars) for v in stars)+sum(tuple(-x for x in v) not in removed for v in removed),'all_squared_norms':32,'dense_minimum_distance':distance,'ordinary_pair_product_upper_bounds':maxima,'certificate_sha256':hashlib.sha256(raw).hexdigest(),'arithmetic':'exact Python integers, individual signed points, no absolute-product test','seconds':time.monotonic()-start}

if __name__=='__main__':
 if not __debug__:raise SystemExit('Run without -O.')
 path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('certificate33_370408_explicit.json')
 print(json.dumps(verify(path),indent=2))
