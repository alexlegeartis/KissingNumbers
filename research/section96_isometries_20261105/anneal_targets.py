"""Minimize the exact number of old-section conflicts over mixed-parity targets."""
import numpy as np,subprocess,re,sys,math
seed=int(sys.argv[1]) if len(sys.argv)>1 else 920201;rng=np.random.default_rng(seed);W=np.array([x for x in range(16) if x.bit_count()==2]);O=np.array([x for x in range(16) if x.bit_count()%2]);target=np.r_[rng.choice(W,4),rng.choice(O,20)].astype(int)
def score(x):
 word=''.join(format(int(v),'x') for v in x);p=subprocess.run(['./exact_count_decoder',word],capture_output=True,text=True)
 if p.returncode:raise RuntimeError(p.stderr)
 if p.stdout.startswith('EXACT_HOLE'):return 0,word
 q=re.search(r'CONFLICTS (\d+)',p.stdout);assert q;return int(q[1]),word
cur,word=score(target);best=cur;print('start',seed,cur,word,flush=True)
for i in range(5000):
 q=target.copy();k=1 if rng.random()<.85 else int(rng.integers(2,5));idx=rng.choice(24,k,replace=False)
 for j in idx:q[j]=rng.choice(W if j<4 else O)
 val,word=score(q);T=2+12*(1-(i%500)/500)
 if val<=cur or rng.random()<math.exp((cur-val)/T):target=q;cur=val
 if val<best:
  best=val;print('best',best,'iteration',i,'word',word,flush=True);open(f'search92_count_anneal_best_{seed}.txt','w').write(word+'\n')
 if val==0:open(f'search92_count_anneal_HOLE_{seed}.txt','w').write(word+'\n');break
 if i%100==99:print('progress',i+1,'current',cur,'best',best,flush=True)
print('DONE best',best,flush=True)
