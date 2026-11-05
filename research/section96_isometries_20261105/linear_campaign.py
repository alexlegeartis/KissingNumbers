"""Schedule distinct 2D parity spaces, retaining individual logs."""
import numpy as np,subprocess,concurrent.futures,sys,pathlib,os
rng=np.random.default_rng(922119);seen=set();pairs=[]
while len(pairs)<int(os.getenv('LINEAR_SPACES','12')):
 u=int(rng.integers(1,256));v=int(rng.integers(1,256))
 if u==v:continue
 key=tuple(sorted((u,v,u^v)))
 if key in seen:continue
 seen.add(key);pairs.append((key[0],key[1]))
def run(item):
 i,(u,v)=item;log=pathlib.Path(f'sectionlinear_campaign_{i:02}_{u}_{v}.log')
 with log.open('w') as f:
  p=subprocess.run(['python','-u','search_linear.py','4','4',str(u),str(v),os.getenv('LINEAR_ROUNDS','60'),str(193000+i)],stdout=f,stderr=subprocess.STDOUT)
 print('space',i,u,v,'exit',p.returncode,flush=True)
 if p.returncode:raise RuntimeError(f'campaign subprocess failed: {log}')
with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.getenv('LINEAR_WORKERS','3'))) as ex:
 list(ex.map(run,enumerate(pairs)))
