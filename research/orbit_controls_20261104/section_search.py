"""Search a common signed-pair normal catalogue, then certify best sections."""
from pathlib import Path
from itertools import combinations
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
from discovery.engines.geometry import section
from discovery.exact_gram import exact_integer_rank

def main():
 out=Path(__file__).resolve().parent;start=time.monotonic()
 source=Artifact.from_dict(json.loads((out/'leech_artifact.json').read_text()))
 context=RunContext(ROOT,{'max_materialized_points':250000,'max_orbit_points':250000,'max_rank_work':1000000000,'max_pairs':250000})
 if verify_artifact(source,context).status!='exact':raise ValueError('Fresh complete Leech replay failed')
 X=np.array(source.data['coordinates'],np.int16);active=np.ones(len(X),bool);normals=[];records=[]
 for depth in range(1,5):
  current=X[active];candidates=[]
  for i,j in combinations(range(24),2):
   for sign in (-1,1):
    normal=[0]*24;normal[i]=1;normal[j]=sign
    rank,_=exact_integer_rank(normals+[normal],context)
    if rank!=depth:continue
    count=int(np.count_nonzero(current[:,i]+sign*current[:,j]==0))
    candidates.append((count,i,j,sign,normal))
  best=max(candidates,key=lambda item:(item[0],-item[1],-item[2],-item[3]))
  count,i,j,sign,normal=best;normals.append(normal);active&=X[:,i]+sign*X[:,j]==0
  result=section([source],{'normals':normals},context)
  if result.facts['cardinality']!=count:raise ValueError('Independent section differs from searched count')
  # Each artifact is self-contained and independently replays the full source.
  path=out/f'section{24-depth}_artifact.json'
  path.write_text(json.dumps(result.to_dict(),separators=(',',':'))+'\n')
  record={'dimension':24-depth,'facts':result.facts,'normals':[x[:] for x in normals],
          'normal_catalogue':'all signed coordinate-pair normals, independent of targets',
          'candidate_normals_tested':len(candidates),'artifact':path.name,'verification':result.verification.to_dict(),
          'new_sota':False}
  records.append(record);print(json.dumps({'dimension':24-depth,'points':count,'tested':len(candidates)}),flush=True)
  (out/'section_search.json').write_text(json.dumps({'scope':'Finite central-section construction search, followed by exact subset proof; no global optimum or new record claim','records':records,'seconds':time.monotonic()-start},indent=2)+'\n')
if __name__=='__main__':main()
