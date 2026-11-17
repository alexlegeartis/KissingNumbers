import numpy as np,json
from pathlib import Path
import sympy as sp
P=Path(__file__).resolve().parent;x=sp.symbols('x');p=2*x**27+2*x**25+x**23+x**21+2*x**20+1;g=sp.Poly(sp.gcd(p,x**28+1,modulus=3),x,modulus=3);coef=np.array([int(g.nth(i))%3 for i in range(g.degree()+1)]);G=np.zeros((20,28),dtype=np.int64)
assert g.degree()==8
for i in range(20):G[i,i:i+9]=coef
R=G[:,:27][:,15:27].T.copy()%3;piv=[];rank=0
for j in range(20):
 ix=np.flatnonzero(R[rank:,j]);
 if not len(ix):continue
 k=rank+ix[0];R[[rank,k]]=R[[k,rank]];R[rank]=R[rank]*int(R[rank,j])%3;c=R[:,j].copy();c[rank]=0;R=(R-c[:,None]*R[rank])%3;piv.append(j);rank+=1
 if rank==12:break
free=[j for j in range(20)if j not in piv];assert len(free)==8;N=np.zeros((8,20),int)
for i,j in enumerate(free):N[i,j]=1;N[i,piv]=-R[:,j]%3
H=N@G%3;assert np.all(H[:,15:27]==0);H=H[:,:15];K=(np.arange(3**8)[:,None]//3**np.arange(8)%3)@H%3;assert len(np.unique(K,axis=0))==6561;wt=np.count_nonzero(K,axis=1);assert wt[0]==0 and min(wt[1:])==5
np.save(P/'phase_generator.npy',H);np.save(P/'phase_code.npy',K.astype(np.uint8));source=Path('research/pless48_fiber9_cert_20261117');outer=np.load(source/'min15.npy');GG=np.load(source/'generator.npy');assert np.all(GG@GG.T%3==0);assert np.all(outer[:,24:]==outer[:,:24]@GG[:,24:]%3);assert np.all(np.count_nonzero(outer,axis=1)==15);assert len(np.unique(outer,axis=0))==415104
profiles=[]
for cs in [list(range(3)),list(range(4)),list(range(5)),list(range(6)),[12,14,25,41,42,44],list(range(7))]:
 sub=outer[np.all(outer[:,cs]==0,axis=1)];body=[i for i in range(48)if i not in cs];sub=sub[:,body];name='outer_'+str(len(body))+'_'+('_'.join(map(str,cs)));np.savez_compressed(P/(name+'.npz'),words=sub.astype(np.uint8));profiles.append({'shortened_coordinates':cs,'body_blocks':len(body),'min15_words':len(sub),'phase_words':6561,'sparse_points':len(sub)*6561,'file':name+'.npz'})
json.dump({'phase_source':'https://www.codetables.de/BKLC/BKLC.php?q=3&n=15&k=8','phase_parent_polynomial':str(p),'phase_gcd_generator_low_to_high':coef.tolist(),'phase_generator':H.tolist(),'phase_rank':8,'phase_size':6561,'phase_min_distance':5,'phase_weight_distribution':dict(zip(*[a.tolist()for a in np.unique(wt,return_counts=True)])),'outer_generator':GG.tolist(),'profiles':profiles,'outer_minimum_premise':'The ternary Pless[48,24,15] code hasminimum15, established in parent full enumeration/source. Distinct selectedweight15 words verified directly.'},open(P/'certificate.json','w'),indent=2);print(json.load(open(P/'certificate.json'))['profiles'])
