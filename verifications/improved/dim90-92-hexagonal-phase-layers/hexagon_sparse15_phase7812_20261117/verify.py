import numpy as np,json,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent;raw=(P/'phase7812.txt').read_bytes();assert hashlib.sha256(raw).hexdigest()=='d5ba3940239da9750271082790a2ac3f62aa52a14c4ee66727f9a66121f593e2';lines=raw.decode().splitlines();assert len(lines)==7812 and all(len(x)==15 and set(x)<=set('012')for x in lines);C=np.array([[int(v)for v in row]for row in lines],dtype=np.uint16);assert len(np.unique(C,axis=0))==7812;one=((C==1)*(1<<np.arange(15))).sum(1).astype(np.uint16);two=((C==2)*(1<<np.arange(15))).sum(1).astype(np.uint16);pc=np.array([i.bit_count()for i in range(1<<15)],dtype=np.uint8);hist=np.zeros(16,np.int64)
for i in range(7811):hist+=np.bincount(pc[(one[i]^one[i+1:])|(two[i]^two[i+1:])],minlength=16)
assert hist.sum()==30509766 and not np.any(hist[:5])and hist[5]>0
c=json.load(open(P/'outer_certificate.json'));G=np.array(c['outer_generator'],dtype=np.int64);assert np.array_equal(G[:,:24],np.eye(24,dtype=int));reports=[]
for r in c['profiles']:
 W=np.load(P/r['file'])['words'].astype(np.int64);assert len(W)==r['min15_words']and len(np.unique(W,axis=0))==len(W)and np.all(np.count_nonzero(W,axis=1)==15);body=[j for j in range(48)if j not in r['shortened_coordinates']];full=np.zeros((len(W),48),np.int64);full[:,body]=W;assert np.array_equal(full[:,24:],full[:,:24]@G[:,24:]%3);assert set(map(bytes,W.astype(np.uint8)))==set(map(bytes,((3-W)%3).astype(np.uint8)));reports.append({'body_blocks':len(body),'outer_words':len(W),'shortened_coordinates':r['shortened_coordinates'],'sparse_points':len(W)*7812})
h=np.array([[3,1],[-3,1],[0,-2]]);t=np.array([[2,0],[-1,1],[-1,-1]]);M=np.diag([1,3]);assert np.all(np.diag(h@M@h.T)==12)and set((h@M@t.T).ravel())=={-6,0,6}
for i in range(3):
 for j in range(3):
  v=h[i]+h[j];assert v@M@v>=12
  if i!=j:v=h[i]-h[j];assert v@M@v==36
base=2*3**21+130944*7812;assert base==21943640934
json.dump({'status':'PASS','phase_code_size':7812,'all_unordered_pairs':int(hist.sum()),'phase_minimum_distance':5,'distance_histogram':hist.tolist(),'outer_profiles':reports,'geometry_exact':True,'dense90_plus_sparse':base,'with_orthogonal_complement91':base+2,'with_orthogonal_complement92':base+6,'outer_minimum15_premise':c['outer_minimum_premise']},open(P/'replay.json','w'),indent=2);print('PASS7812phaseallpairs,outermembership/counts,geometry;90/91/92',base,base+2,base+6)
