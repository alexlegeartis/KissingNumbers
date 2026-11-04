"""Both sign parities, partial stars, shared root deletions via the engine."""
from pathlib import Path
from itertools import combinations,product
import argparse,json,subprocess,sys,tempfile,time,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.characters import assemble_orbits
PACKAGE=ROOT/'verifications/improved/dim32-35-signed-layers'
OUT=Path(__file__).resolve().parent

def key(v):return min(tuple(v),tuple(-x for x in v))
def star_points(n,a,leaves,parity=None):
 for signs in product((-1,1),repeat=4):
  if parity is not None and sum(x<0 for x in signs)%2!=parity:continue
  v=[0]*n;v[a]=4
  for b,s in zip(leaves,signs):v[b]=2*s
  yield v

def search(n,seconds):
 start=time.monotonic()
 if n==39:
  source=ROOT/'verifications/improved/dim39-deformed-orbits/source'
  basepath=source/'orbit_certificate.json';dense=json.loads(basepath.read_text());sparse=json.loads((source/'base_certificate.json').read_text())
  ss=[tuple(j for j,x in enumerate(line) if x=='1') for line in (source/'source/supports_3424.txt').read_text().splitlines()]
  for change in sparse['support_replacements']:
   if list(ss[change['remove_index']])!=change['remove']:raise ValueError('Sparse replacement source')
   ss[change['remove_index']]=tuple(change['add'])
  masks=[sum(1<<j for j in s) for s in ss]
  c={'dimension':39,'full_even_supports':masks,'points':len(dense['representatives'])*len(dense['group_words'])+128*len(masks)+2*n*(n-1)}
  removed_pairs={tuple(sorted((s['center'],b))) for s in sparse['heavy_stars'] for b in s['leaves']}
  old={'stars':sparse['heavy_stars'],'removed_root_lines':[list(p)+( [s]) for p in sorted(removed_pairs) for s in (-1,1)]}
 else:
  basepath=PACKAGE/'source'/('multi35/certificate35.json' if n==35 else f'classic/certificate{n}.json')
  c=json.loads(basepath.read_text());old=json.loads((PACKAGE/f'root_trades/certificate{n}.json').read_text())
 exceptional=np.array([[(-2 if s>>j&1 else 2) if m>>j&1 else 0 for j in range(n)] for m,s in c.get('exceptional_antipodal_lines',[])],np.int64).reshape((-1,n))
 heads=[];removed=set()
 if n==34:
  a,b=c['heavy_pair_extension']['coordinates'];removed.add((a,b,-1))
  for t in (-1,1):
   v=[1]*n;v[a],v[b]=3*t,-3*t;heads.append(v)
 elif n==35:heads=c['head_antipodal_lines'];removed=set(map(tuple,c['removed_root_lines']))
 heads=np.array(heads,np.int64).reshape((-1,n))
 supports=c.get('full_even_supports',[m for m,p in c.get('full_parity_supports',[])])
 with tempfile.TemporaryDirectory(prefix='partial-star-') as temp:
  temp=Path(temp);(temp/'supports').write_text(''.join(f'{s}\n' for s in supports))
  subprocess.run(['g++','-O3','-std=c++17',str(PACKAGE/'root_trades/enumerate_stars.cpp'),'-o',str(temp/'enum')],check=True)
  subprocess.run([str(temp/'enum'),str(n),str(temp/'supports'),str(temp/'stars')],check=True)
  stars=[tuple(map(int,line.split())) for line in (temp/'stars').read_text().splitlines()]
 old_roots={}
 for a,b in combinations(range(n),2):
  for sign in (-1,1):
   if (a,b,sign) in removed:continue
   v=[0]*n;v[a]=4;v[b]=4*sign;old_roots[(a,b,sign)]=key(v)
 candidates={};origins={}
 for a,mask in stars:
  leaves=[j for j in range(n) if mask>>j&1]
  for v in star_points(n,a,leaves):
   V=np.array(v,np.int64)
   if len(exceptional) and np.max(np.abs(exceptional@V))>16:continue
   if len(heads) and np.any(4*(heads@V)**2>32*np.sum(heads*heads,axis=1)):continue
   k=key(v);candidates[k]=v;origins[k]={'center':a,'leaves':leaves,'leaf_parity':sum(v[b]<0 for b in leaves)%2}
 initial={v for line,v in old_roots.items() if line not in set(map(tuple,old['removed_root_lines']))}
 for star in old['stars']:
  for v in star_points(n,star['center'],star['leaves'],0):
   k=key(v)
   if k not in candidates:raise ValueError('Current star point absent from candidate grammar')
   initial.add(k)
 shapes=sorted(set(old_roots.values())|set(candidates));index={v:i for i,v in enumerate(shapes)}
 params={'generators':[[1]*n],'shapes':[list(v) for v in shapes],'initial_shapes':[index[v] for v in sorted(initial)],'seconds':seconds}
 context=RunContext(ROOT,{'max_gram_entries':2*len(shapes)**2+1,'max_rank_work':len(shapes)**2*(n+3)+1,'max_pairs':100000000,'max_materialized_points':max(20000,len(shapes)*n),'solver_seconds':max(30,seconds)})
 print(json.dumps({'dimension':n,'support_stars':len(stars),'candidate_lines':len(candidates),'root_lines':len(old_roots),'initial_component_points':2*len(initial)}),flush=True)
 artifact=assemble_orbits([],params,context)
 chosen={tuple(row) for row in artifact.data['coordinates']};new={key(row) for row in artifact.data['coordinates'] if key(row) in candidates};kept={line for line,v in old_roots.items() if v in chosen or tuple(-x for x in v) in chosen}
 fixed=c['points']-2*len(old_roots);total=fixed+artifact.facts['cardinality'];old_total=fixed+2*len(initial)
 report={'dimension':n,'strongest_captured_target':{33:370376,34:390990,35:415334,39:769180}[n],'baseline_points':c['points'],'retained_fixed_points':fixed,'initial_total':old_total,'proposed_total':total,'gain':total-old_total,'star_supports':len(stars),'candidate_antipodal_lines':len(candidates),'selected_new_lines':[list(v) for v in sorted(new)],'removed_root_lines':[list(line) for line in sorted(set(old_roots)-kept)],'search':artifact.data['search'],'facts':artifact.facts,'budget':context.budget,'source_baseline':str(basepath.relative_to(ROOT)),'source_baseline_sha256':hashlib.sha256(basepath.read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'new_sota':False,'scope':'Only the root/new-point component is exact here. Fixed dense/middle/heads require independent full composition proof before any lower-bound claim.'}
 (OUT/f'partial_stars{n}.json').write_text(json.dumps(report,indent=2)+'\n')
 (OUT/f'partial_stars{n}_artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 print(json.dumps({k:report[k] for k in ('dimension','initial_total','proposed_total','gain','seconds')}),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('dimensions',nargs='*',type=int,default=[33,34,35]);ap.add_argument('--seconds',type=float,default=30);args=ap.parse_args()
 for n in args.dimensions:search(n,args.seconds)
