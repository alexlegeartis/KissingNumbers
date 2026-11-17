import numpy as np,json
from pathlib import Path
P=Path(__file__).resolve().parent;c=json.load(open(P/'certificate.json'));G=np.array(c['outer_generator'],dtype=np.int64);H=np.array(c['phase_generator'],dtype=np.int64);K=(np.arange(6561)[:,None]//3**np.arange(8)%3)@H%3;assert len(np.unique(K,axis=0))==6561 and min(np.count_nonzero(K[1:],axis=1))==5
# Coordinates(a,b) mean(a/2,b*sqrt3/2), so four times inner product is integral.
t=np.array([[2,0],[-1,1],[-1,-1]]);h=np.array([[3,1],[-3,1],[0,-2]]);M=np.diag([1,3]);assert np.all(np.diag(h@M@h.T)==12);assert set((h@M@t.T).ravel())=={-6,0,6};assert set((h@M@h.T).ravel())=={-6,12}
for i in range(3):
 for j in range(3):
  d=h[i]+h[j];assert int(d@M@d)>=12
  if i!=j:d=h[i]-h[j];assert int(d@M@d)==36
reports=[]
for row in c['profiles']:
 W=np.load(P/row['file'])['words'].astype(np.int64);assert len(W)==row['min15_words']and len(np.unique(W,axis=0))==len(W)and np.all(np.count_nonzero(W,axis=1)==15);full=np.zeros((len(W),48),np.int64);body=[i for i in range(48)if i not in row['shortened_coordinates']];full[:,body]=W;assert np.array_equal(full[:,24:],full[:,:24]@G[:,24:]%3);assert len(W)*len(K)==row['sparse_points'];neg=(3-W)%3;assert set(map(bytes,W.astype(np.uint8)))==set(map(bytes,neg.astype(np.uint8)));reports.append(row)
json.dump({'status':'PASS','phase_words':6561,'phase_minimum':5,'coordinate_distance_tables_exact':True,'profiles':reports,'parent_minimum15_premise':c['outer_minimum_premise'],'dense90_plus_sparse45':2*3**21+130944*6561},open(P/'verification.json','w'),indent=2);print('PASS phasecode,geometrytables,outerwordmembership/counts,antipodality')
