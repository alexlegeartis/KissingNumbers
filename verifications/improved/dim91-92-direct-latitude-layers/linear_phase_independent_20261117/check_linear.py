from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent;S=P/'linear_phases'
reports={}
for n,k in [(14,8),(13,7)]:
 G=np.load(S/f'generator{n}.npy').astype(np.int64);assert G.shape==(k,n)
 expected=np.zeros((k,n),dtype=np.int64)
 for i in range(k):expected[i,i:i+7]=[1,1,0,1,0,1,1]
 assert np.array_equal(G,expected)
 # Generate the vector space by iterative closure, independently of base3 message multiplication.
 words=np.zeros((1,n),dtype=np.int64)
 for row in G:words=np.concatenate((words,(words+row)%3,(words+2*row)%3))
 assert len(set(map(bytes,words.astype(np.uint8))))==3**k
 weights=np.count_nonzero(words,axis=1);hist=np.bincount(weights,minlength=n+1)
 assert hist[0]==1 and hist[1:5].sum()==0 and hist[5]>0
 saved=np.load(S/f'phase{n}.npy');assert set(map(bytes,saved))==set(map(bytes,words.astype(np.uint8)))
 reports[str(n)]=dict(dimension=k,size=len(words),minimum=5,histogram=hist.tolist())
reports['total91']=21943640934+2*31680*6561+2
reports['total92']=21943640934+4*31680*6561+4*6930*2187+6
assert reports['total91']==22359345896 and reports['total92']==22835674500
(P/'linear_result.json').write_text(json.dumps(reports,indent=2));print(json.dumps(reports,indent=2))
