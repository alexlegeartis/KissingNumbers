"""Allow individual signed stars and individual root deletions above 31D."""
import argparse,itertools,json,subprocess,sys,time
from pathlib import Path
import numpy as np
sys.path.append('/workspace/work/deps')
from ortools.sat.python import cp_model
ROOT=Path(__file__).resolve().parents[2]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=180);ap.add_argument('--seed',type=int,default=811);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic()
 pkg=ROOT/'verifications/improved/dim32-35-signed-layers';c=json.loads((pkg/'source/multi35/certificate35.json').read_text());old=json.loads((pkg/'root_trades/certificate35.json').read_text());d=35
 support=a.out/'supports.txt';support.write_text(''.join(f'{mask}\n' for mask,parity in c['full_parity_supports']))
 binary=a.out/'enumerate';subprocess.run(['c++','-O3','-std=c++17',str(pkg/'root_trades/enumerate_stars.cpp'),'-o',str(binary)],check=True);subprocess.run([str(binary),str(d),str(support),str(a.out/'stars.txt')],check=True)
 heads=np.array(c['head_antipodal_lines'],np.int64);hn=np.sum(heads*heads,axis=1);V=[];owner=[]
 for line in (a.out/'stars.txt').read_text().splitlines():
  center,mask=map(int,line.split());leaves=[i for i in range(d) if mask>>i&1]
  for signs in itertools.product((-1,1),repeat=5):
   v=np.zeros(d,np.int16);v[center]=4*signs[0];v[leaves]=2*np.array(signs[1:]);dot=v.astype(np.int64)@heads.T
   if np.any(dot*dot>8*hn):continue
   V.append(v);owner.append((center,leaves))
 V=np.unique(np.array(V,np.int16),axis=0);lookup={v.tobytes():i for i,v in enumerate(V)}
 initial=set()
 for star in old['stars']:
  for cs in (-1,1):
   for ss in itertools.product((-2,2),repeat=4):
    if sum(x<0 for x in ss)%2:continue
    v=np.zeros(d,np.int16);v[star['center']]=4*cs;v[star['leaves']]=ss;initial.add(lookup[v.tobytes()])
 fixed_removed=set(map(tuple,c['removed_root_lines']));blockers=[];roots={}
 for v in V:
  center=int(np.flatnonzero(abs(v)==4)[0]);bs=[]
  for leaf in np.flatnonzero(abs(v)==2):
   i,j=sorted((center,int(leaf)));si,sj=int(np.sign(v[i])),int(np.sign(v[j]))
   if (i,j,si*sj) in fixed_removed:continue
   key=(i,j,si,sj);roots.setdefault(key,len(roots));bs.append(roots[key])
  blockers.append(bs)
 original_removed=set(b for i in initial for b in blockers[i]);assert len(initial)-len(original_removed)==old['points']-c['points']==20
 model=cp_model.CpModel();x=[model.new_bool_var(f'x{i}') for i in range(len(V))];y=[model.new_bool_var(f'r{i}') for i in roots];edges=0
 for i,bs in enumerate(blockers):
  for b in bs:model.add(x[i]<=y[b])
 F=V.astype(np.float32)
 for at in range(0,len(V),96):
  bad=F[at:at+96]@F.T>16
  for i,j in zip(*np.nonzero(bad)):
   i+=at
   if i<j:model.add_at_most_one(x[int(i)],x[int(j)]);edges+=1
 score=sum(x)-sum(y);model.add(score>=21);model.maximize(score)
 for i,v in enumerate(x):model.add_hint(v,i in initial)
 for b,v in enumerate(y):model.add_hint(v,b in original_removed)
 print('MODEL35',len(V),'roots',len(roots),'edges',edges,'seconds',time.monotonic()-start,flush=True)
 solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.random_seed=a.seed;solver.parameters.log_search_progress=True
 st=solver.solve(model);report=dict(status=solver.status_name(st),dimension=d,candidates=len(V),edges=edges,baseline=old['points'],seconds=time.monotonic()-start)
 if st in [cp_model.FEASIBLE,cp_model.OPTIMAL]:
  take=np.flatnonzero([solver.boolean_value(v) for v in x]);removed={b for i in take for b in blockers[i]};S=V[take].astype(np.int64);G=S@S.T;np.fill_diagonal(G,-32);assert G.max()<=16
  rootlist=[k for k,i in roots.items() if i in removed];count=c['points']+len(take)-len(removed);assert count>old['points']
  proposal=dict(points=count,baseline=str(pkg/'source/multi35/certificate35.json'),signed_star_points=S.tolist(),removed_root_points=rootlist)
  (a.out/f'proposal35_{count}.json').write_text(json.dumps(proposal,indent=2)+'\n');report.update(count=count,gain=count-old['points'])
 (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
