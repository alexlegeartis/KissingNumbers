"""Search an index-eight DGLS subcode under independent block isometries."""
import argparse,itertools,json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.block_isometry_search import search_block_isometries
p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=60);p.add_argument('--rounds',type=int,default=8);p.add_argument('--output',type=pathlib.Path);args=p.parse_args()
g=[7,1,4,4,2,7,3,3,4,7,2,1];G=[[0]*24 for _ in range(12)]
for i in range(12):
 for j,x in enumerate(g):G[i][i+j]=x
 G[i][23]=3
assert all(G[i][j]==0 for i in range(12) for j in range(i))
assert all(G[i][i]%2==1 for i in range(12))
base=[row[:] for row in G[:3]]
for i in range(3):G[i]=[2*x%8 for x in G[i]]
offsets=[[sum(((r>>i)&1)*base[i][j] for i in range(3))%8 for j in range(24)] for r in range(8)]
gray=[sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)]
pairs=[(0,15),(3,12),(5,10),(6,9)];menu=[]
for perm in itertools.permutations(range(4)):
 for flips in range(16):
  mp={}
  for i,(a,b) in enumerate(pairs):
   v=pairs[perm[i]];f=(flips>>i)&1;mp[a]=v[f];mp[b]=v[1-f]
  menu.append([mp[x] for x in gray])
assert len({tuple(x) for x in menu})==384
params={'generator':G,'modulus':8,'bit_width':4,'block_map':gray,'map_menu':menu,'groups':list(range(24)),'references':[[0]*24 for _ in range(8)],'reference_offsets':offsets,'minimum_distance':24,'rounds':args.rounds,'oracle_steps':500,'exact_pair_limit':4096,'seed':9618,'seconds':args.seconds}
ctx=RunContext(ROOT,{'max_trials':max(200,args.rounds),'solver_seconds':max(args.seconds,1),'max_rank_work':1000000000})
r=search_block_isometries(params,ctx);print(json.dumps({'status':r['status'],'oracle':r['oracle'],'witnesses':len(r['rejection_witnesses']),'unresolved_candidates':len(r['candidates']),'scope':r['scope']},indent=2))
if args.output:args.output.write_text(json.dumps(r))
