from pathlib import Path
import numpy as np,json,hashlib
P=Path(__file__).resolve().parent
C=json.loads((P/'certificate.json').read_text());G=np.load(P/'generator.npy').astype(np.int64)
assert G.shape==(24,48) and np.array_equal(G[:,:24],np.eye(24,dtype=int))
K=np.load(P/'source_phase7812.npy');H=np.load(P/'phase14.npy')
assert np.array_equal(H,np.delete(K[K[:,14]==1],14,axis=1))
assert H.shape==(2609,14) and np.all(H<3) and len(np.unique(H,axis=0))==2609
hist=np.zeros(15,np.int64)
for i in range(len(H)-1):hist+=np.bincount(np.count_nonzero(H[i+1:]!=H[i],axis=1),minlength=15)
assert hist[:5].sum()==0 and hist[5]>0 and hist.sum()==2609*2608//2
assert hist.tolist()==C['phase_distance_histogram']
for row in C['outer_layers']:
 B=np.load(P/row['file']);assert B.shape==(31680,45) and np.all(B<3)
 assert np.all(np.count_nonzero(B,axis=1)==14) and len(np.unique(B,axis=0))==31680
 F=np.column_stack((np.tile(row['prefix'],(len(B),1)),B)).astype(np.int64)
 assert np.array_equal((F[:,:24]@G)%3,F)
# Coordinates encode (a/2,b*sqrt(3)/2), so dot numerator is aa'+3bb'.
R=np.array([[3,1],[-3,1],[0,-2]],dtype=np.int64)
T=np.array([[2,0],[-1,1],[-1,-1]],dtype=np.int64)
def dot(A,B):return A[:,0,None]*B[None,:,0]+3*A[:,1,None]*B[None,:,1]
assert np.array_equal(dot(R,R),np.array([[12,-6,-6],[-6,12,-6],[-6,-6,12]]))
assert np.max(np.abs(dot(R,T)))==6
assert np.min(24+2*dot(R,R))==12 # opposite signs, four times squared distance
assert len(np.unique(np.vstack((R,-R)),axis=0))==6
assert 3*15==45 and 9*5==45 and 3*13+6==45
assert 45+42-3*14==45 and 4*135<2025
assert C['total91']==21943640934+2*31680*2609+2==22108947176
assert C['total92']==21943640934+4*31680*2609+6==22274253420
print('PASS: finite inputs, exact local geometry, counts. Pless minimum15 and equator validity remain stated premises.')
