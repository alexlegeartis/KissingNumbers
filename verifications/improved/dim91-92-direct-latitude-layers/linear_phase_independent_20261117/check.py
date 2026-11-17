from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent;S=P/'input14';G=np.load(S/'generator.npy').astype(np.int64)
assert G.shape==(24,48) and np.array_equal(G[:,:24],np.eye(24,dtype=np.int64)) and not np.any(G@G.T%3)
assert hashlib.sha256((P/'prior_sparse15_certificate.json').read_bytes()).hexdigest()=='5bddcfc35bb30cec8257e7b1a30775f3153fc1b796a9d12cb087b34265398eb4'
oldcert=json.load(open(P/'prior_sparse15_certificate.json'));old=np.array(oldcert['outer_generator'],dtype=np.int64);assert np.array_equal(old,G)
H=np.load(S/'phase14.npy');K=np.load(S/'source_phase7812.npy');assert H.shape==(2609,14) and np.all((H>=0)&(H<=2));assert np.array_equal(H,np.delete(K[K[:,14]==1],14,axis=1))
# Independently count all pairs by packed two-bit symbol XOR, not array Hamming reductions.
words=[]
for row in H:words.append(sum(int(v)<<(2*j) for j,v in enumerate(row)))
assert len(set(words))==2609;mask=sum(1<<(2*j)for j in range(14));hist=[0]*15
for i,a in enumerate(words):
 for b in words[i+1:]:
  x=a^b;hist[((x|(x>>1))&mask).bit_count()]+=1
assert sum(hist)==3402136 and sum(hist[:5])==0
C=json.load(open(S/'certificate.json'));assert hist==C['phase_distance_histogram']
for row in C['outer_layers']:
 B=np.load(S/row['file']).astype(np.int64);assert B.shape==(31680,45) and np.all((B>=0)&(B<=2));assert np.all(np.count_nonzero(B,axis=1)==14)
 assert len({bytes(r) for r in B.astype(np.uint8)})==31680
 full=np.concatenate([np.broadcast_to(row['prefix'],(len(B),3)),B],axis=1)
 assert np.array_equal(full[:,:24]@G%3,full)
# Dot products in coordinates (a/2,b*sqrt(3)/2), evaluated as integers/4.
R=[(3,1),(-3,1),(0,-2)];T=[(2,0),(-1,1),(-1,-1)]
dot=lambda a,b:a[0]*b[0]+3*a[1]*b[1]
assert all(dot(x,x)==12 for x in R)
assert all(dot(R[i],R[j])==-6 for i in range(3) for j in range(i))
assert max(abs(dot(x,t)) for x in R for t in T)==6
assert min(dot((a[0]+b[0],a[1]+b[1]),(a[0]+b[0],a[1]+b[1]))for a in R for b in R)==12
assert 135*4<2025
layers=31680*2609;equator=21943640934
report=dict(status='FINITE INPUTS AND UNION GEOMETRY PASS',phase_pairs=sum(hist),phase_minimum=next(i for i,a in enumerate(hist)if a),each_outer=31680,points_per_layer=layers,total91=equator+2*layers+2,total92=equator+4*layers+6,source_generator_sha256=hashlib.sha256((S/'generator.npy').read_bytes()).hexdigest(),premises=['Previously established Pless generator minimum15','Previously audited 21943640934-point equatorial union'])
assert report['total91']==22108947176 and report['total92']==22274253420
(P/'result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

# Independent full finite checks of the four weight2-prefix corner layers.
import itertools
S2=P/'input13';assert np.array_equal(np.load(S2/'generator.npy'),G)
J=np.load(S2/'phase13.npy');KK=np.load(S2/'source_phase7812.npy')
assert np.array_equal(KK,K) and np.array_equal(J,np.delete(KK[np.all(KK[:,[0,3]]==[0,2],axis=1)],[0,3],axis=1))
assert J.shape==(888,13) and np.all((J>=0)&(J<=2))
codes=[sum(int(v)<<(2*j)for j,v in enumerate(row))for row in J];assert len(set(codes))==888
mask13=sum(1<<(2*j)for j in range(13));hist13=[0]*14
for i,a in enumerate(codes):
 for b in codes[i+1:]:
  t=a^b;hist13[((t|(t>>1))&mask13).bit_count()]+=1
assert sum(hist13)==888*887//2 and not any(hist13[:5])
for x,y in itertools.product([1,2],repeat=2):
 B=np.load(S2/f'outer_{x}{y}0.npy').astype(np.int64);assert B.shape==(6930,45) and np.all((B>=0)&(B<=2))
 assert np.all(np.count_nonzero(B,axis=1)==13) and len(set(map(bytes,B.astype(np.uint8))))==6930
 full=np.concatenate([np.broadcast_to([x,y,0],(6930,3)),B],axis=1);assert np.array_equal(full[:,:24]@G%3,full)
f=[0,1,-1]
for p in itertools.product(range(3),repeat=2):
 for q in itertools.product(range(3),repeat=2):
  assert sum((f[a]-f[b])**2for a,b in zip(p,q))>=sum(a!=b for a,b in zip(p,q))
assert 4*270<2025
report['corner_phase_pairs']=sum(hist13);report['corner_added']=4*6930*888;report['total92_with_corners']=report['total92']+report['corner_added'];assert report['total92_with_corners']==22298868780
(P/'result.json').write_text(json.dumps(report,indent=2));print('GRID EXTENSION PASS',json.dumps(report,indent=2))
