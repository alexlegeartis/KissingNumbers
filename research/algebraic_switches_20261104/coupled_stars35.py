"""Joint middle-family/root-star replacement through the generic selector.

The conservative middle oracle deletes an entire 128-point sign family if
any of its points can conflict. Both leaf parities are searched. This extends
the old fixed-middle star family, not the whole kissing construction space.
"""
from pathlib import Path
from itertools import combinations,product
import argparse,gzip,hashlib,json,subprocess,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.runner import run_program,output_artifact

def vectors(n,center,leaves,parity):
 result=[]
 for s in (-1,1):
  for word in range(16):
   if word.bit_count()%2!=parity:continue
   row=np.zeros(n,np.int64);row[center]=4*s
   row[leaves]=[-2 if word>>k&1 else 2 for k in range(4)];result.append(row)
 return np.array(result,np.int64)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--middle-limit',type=int,choices=(1,2),default=1)
 parser.add_argument('--star-limit',type=int,default=6000);parser.add_argument('--out',type=Path);args=parser.parse_args()
 if args.star_limit<100:raise ValueError('Star budget too small to retain the zero-blocker frontier')
 n=35;here=Path(__file__).resolve().parent;out=args.out or here/('coupled_stars35' if args.middle_limit==1 else 'coupled_stars35_m2');out.mkdir(exist_ok=False);start=time.monotonic()
 package=ROOT/'verifications/improved/dim32-35-signed-layers'
 source=package/'source/multi35/certificate35.json';base=json.loads(source.read_text())
 supports=[int(m) for m,p in base['full_parity_supports']]
 root_removed=set(map(tuple,base['removed_root_lines']));heads=np.array(base['head_antipodal_lines'],np.int64)
 root_lines=[(a,b,s) for a,b in combinations(range(n),2) for s in (-1,1) if (a,b,s) not in root_removed]
 root_index={line:i for i,line in enumerate(root_lines)};middle_offset=len(root_lines)
 costs=[2]*len(root_lines)+[128]*len(supports)
 with tempfile.TemporaryDirectory(prefix='coupled-star-enum-') as temp:
  temp=Path(temp);(temp/'supports').write_text(''.join(str(m)+'\n' for m in supports))
  subprocess.run(['g++','-O3','-std=c++17',str(here/'near_star_supports.cpp'),'-o',str(temp/'enum')],check=True)
  subprocess.run([str(temp/'enum'),str(n),str(temp/'supports'),str(temp/'stars'),str(args.middle_limit)],check=True)
  descriptions=[list(map(int,line.split())) for line in (temp/'stars').read_text().splitlines()]
 stars=[];blockers=[]
 for a,mask,*bad in descriptions:
  leaves=[j for j in range(n) if mask>>j&1]
  for parity in (0,1):
   H=vectors(n,a,leaves,parity)
   if np.any(4*(H@heads.T)**2>32*np.sum(heads*heads,axis=1)):continue
   needed={root_index[line] for b in leaves for s in (-1,1) if (line:=(min(a,b),max(a,b),s)) in root_index}
   needed.update(middle_offset+j for j in bad)
   stars.append({'center':a,'leaf_mask':mask,'leaves':leaves,'parity':parity,'middle_blockers':bad})
   blockers.append(sorted(needed))
 untruncated=len(stars)
 if len(stars)>args.star_limit:
  rng=np.random.default_rng(20261104);order=rng.permutation(len(stars)).tolist();order.sort(key=lambda i:len(stars[i]['middle_blockers']))
  order=order[:args.star_limit];stars=[stars[i] for i in order];blockers=[blockers[i] for i in order]
 print(json.dumps({'support_descriptions':len(descriptions),'admissible_star_families':untruncated,'tested_families':len(stars)}),flush=True)
 conflicts=[]
 for i,A in enumerate(stars):
  a=A['center'];ma=A['leaf_mask']
  for j,B in enumerate(stars[:i]):
   b=B['center'];mb=B['leaf_mask'];overlap=(ma&mb).bit_count()
   forbidden=overlap>0 if a==b else 8*int(bool(mb>>a&1))+8*int(bool(ma>>b&1))+4*overlap>16
   if forbidden:conflicts.append([j,i])
 # Audit the closed pair oracle against independent expanded products before
 # using it; a positive output also requires a complete expanded replay.
 rng=np.random.default_rng(20261104)
 for i,j in rng.integers(0,len(stars),size=(2000,2)):
  if i==j:continue
  A,B=stars[int(i)],stars[int(j)];H=vectors(n,A['center'],A['leaves'],A['parity']);K=vectors(n,B['center'],B['leaves'],B['parity'])
  overlap=(A['leaf_mask']&B['leaf_mask']).bit_count();a,b=A['center'],B['center']
  oracle=overlap>0 if a==b else 8*int(bool(B['leaf_mask']>>a&1))+8*int(bool(A['leaf_mask']>>b&1))+4*overlap>16
  if bool(np.max(abs(H@K.T))>16)!=oracle:raise ValueError('Closed sign-maximum oracle disagrees with expanded points')
 params={'blockers':blockers,'candidate_weights':[16]*len(stars),'removal_costs':costs,
  'candidate_conflicts':conflicts,'seconds':90}
 program={'schema_version':1,'name':'Joint signed-star/middle-family replacement',
  'nodes':[{'id':'exchange','op':'selection.blocker_exchange','params':params}],'outputs':['exchange']}
 receipt=run_program(program,RunContext(ROOT,{'solver_seconds':90,'max_gram_entries':20000000}));proposal=output_artifact(receipt)
 if args.middle_limit==1:
  (out/'program.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
  (out/'receipt.json').write_text(json.dumps(receipt,separators=(',',':'))+'\n')
 else:
  for name,value in [('program.json.gz',program),('receipt.json.gz',receipt)]:
   with gzip.open(out/name,'wt') as f:json.dump(value,f,separators=(',',':'));f.write('\n')
 selected=proposal.data['selected_candidate_indices'];removed=proposal.data['removed_source_indices']
 report={'dimension':n,'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
  'base_points':base['points'],'current_repository':415334,'star_families':len(stars),'star_conflicts':len(conflicts),
  'middle_blocker_limit':args.middle_limit,'untruncated_star_families':untruncated,'star_limit':args.star_limit,
  'star_catalogue_complete':untruncated==len(stars),
  'proposed_total':base['points']+proposal.facts['gain'],'selected_stars':[stars[i] for i in selected],
  'retained_best_total':max(415334,base['points']+proposal.facts['gain']),
  'removed_root_lines':[root_lines[j] for j in removed if j<middle_offset],
  'removed_middle_support_indices':[j-middle_offset for j in removed if j>=middle_offset],
  'solver':proposal.data['solver'],'verification_status':'unverified','new_sota':False,
  'seconds':time.monotonic()-start,'scope':'Complete-family conservative middle deletions, fixed dense code and heads, both leaf parities; finite numerical search, no upper-bound claim. Fresh independent baseline and expanded changed geometry required for any gain.'}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
