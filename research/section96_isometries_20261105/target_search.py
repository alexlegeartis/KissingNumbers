"""Find section targets using exact 4C orbit cuts; proposals are not proofs."""
import numpy as np,subprocess,re,json,os,pathlib
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csr_matrix
cp=pathlib.Path(os.getenv('TARGET92_CHECKPOINT','search92_target_checkpoint.json'));state=json.loads(cp.read_text()) if cp.exists() else {'sources':[],'iteration':0};rng=np.random.default_rng(920021+state['iteration']);g=np.array([7,1,4,4,2,7,3,3,4,7,2,1]);G=np.zeros((12,24),int)
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
Q=np.arange(4096,dtype=np.uint16);bits=np.stack([(Q>>i)&1 for i in range(12)],axis=1).astype(int);golay=(bits@G)%2;assert len(np.unique(golay,axis=0))==4096 and np.min(golay[1:].sum(1))>=8
gray=np.array([sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)]);wt=np.array([int(x).bit_count() for x in range(16)]);menu=[np.array([v for v in range(16) if (wt[v]==2 if j<4 else wt[v]%2==1)],np.uint8) for j in range(24)];offset=np.cumsum([0]+list(map(len,menu)));N=offset[-1];rows=[];lb=[];ub=[]
for j in range(24):
 r=np.zeros(N,np.uint8);r[offset[j]:offset[j+1]]=1;rows.append(r);lb.append(1);ub.append(1)
def add_source(word):
 assert len(word)==12 and set(word)<=set('01234567');b=np.array(list(word),int);source=(b@G)%8;sourceword=gray[source];assert all(wt[sourceword[j]]==2 for j in range(4))
 block=np.zeros((4096,N),np.uint8)
 for j in range(24):
  vals=sourceword[j]^(15*golay[:,j]);block[:,offset[j]:offset[j+1]]=wt[vals[:,None]^menu[j][None,:]]
 rows.extend(block);lb.extend([24]*4096);ub.extend([np.inf]*4096)
for w in state['sources']:add_source(w)
print('binary choices',N,'restored orbits',len(state['sources']),flush=True)
for _ in range(int(os.getenv('TARGET92_ROUNDS','100'))):
 iteration=state['iteration'];res=milp(rng.normal(size=N),integrality=np.ones(N),bounds=Bounds(0,1),constraints=LinearConstraint(csr_matrix(np.asarray(rows)),lb,ub),options={'time_limit':120})
 if res.x is None:print('no incumbent',iteration,res.message,flush=True);break
 rounded=np.rint(res.x).astype(np.int16);assert np.all((rounded==0)|(rounded==1));checks=csr_matrix(np.asarray(rows))@rounded;assert np.all(checks>=np.asarray(lb)-1e-8) and np.all(checks<=np.asarray(ub)+1e-8)
 target=np.array([menu[j][np.argmax(res.x[offset[j]:offset[j+1]])] for j in range(24)],np.uint8);word=''.join(f'{q:x}' for q in target)
 p=subprocess.run(['./exact_target_decoder',word],text=True,capture_output=True)
 if p.returncode:raise RuntimeError((p.returncode,p.stderr))
 match=re.search(r'WITNESS B ([0-7]{12}) distance (\d+)',p.stdout)
 if not match:
  assert p.stdout.startswith('EXACT_HOLE');print('EXACT DECODER HOLE',iteration,word,p.stdout.strip(),'independent replay required',flush=True);pathlib.Path(f'search92_exact_hole_candidate_{iteration}.json').write_text(json.dumps({'target':word,'scope':'exact decoder output; independent replay and literature comparison required'}));break
 else:
  w=match[1];b=np.array(list(w),int);sourceword=gray[(b@G)%8];d=int(wt[sourceword^target].sum());assert d==int(match[2]) and d<24
  
  if len(state['sources'])>=int(os.getenv('TARGET92_MAX_ORBITS','50')):print('orbit cut budget reached',flush=True);break
  add_source(w);state['sources'].append(w);print('iteration',iteration,'bad distance',d,'source orbit cuts',len(state['sources']),flush=True)
 state['iteration']+=1;tmp=cp.with_suffix('.tmp');tmp.write_text(json.dumps(state));tmp.replace(cp)
