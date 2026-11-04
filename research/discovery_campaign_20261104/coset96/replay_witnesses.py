"""Recompute actual rejection pairs; no search imports or sampled-distance claims.
NumPy is used only to read compressed arrays. All witness arithmetic uses Python ints.
"""
from pathlib import Path
import json,time,hashlib
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
start=time.monotonic();source=ROOT/'verifications/improved/dim88-96-z8-golay/data/z8_generator.json';G=json.loads(source.read_text())['generator']
def require(ok,message):
 if not ok:raise ValueError(message)
def vector(v,n,bound,name):
 require(isinstance(v,list) and len(v)==n and all(type(x)is int and 0<=x<bound for x in v),name+' width/domain')
def word(c):
 vector(c,12,8,'message')
 return [sum(a*G[i][j]for i,a in enumerate(c))%8 for j in range(24)]
def gray(a):return sum((((a>>2)&1)^((a&1)*(x&1))^(((a>>1)&1)*((x>>1)&1)))<<x for x in range(4))
def homogeneous(a):return 0 if a==0 else 4 if a==4 else 2
T=[gray(a)for a in range(8)]
for a in range(8):
 for b in range(8):require((T[a]^T[b]).bit_count()==homogeneous((a-b)%8),'Gray metric mismatch')
record=json.loads((OUT/'results.json').read_text())
require(record['generator_sha256']==hashlib.sha256(source.read_bytes()).hexdigest(),'generator source changed')
data=np.load(OUT/'witnesses.npz');ring=0;binary=0
require(len(data['ring_offsets'])==len(data['ring_messages'])==len(data['ring_distances'])==sum(r['candidates']for r in record['ring_results'])+record['ascent_candidates'],'ring witness coverage')
require(len(data['binary_masks'])==len(data['binary_messages'])==len(data['binary_layers'])==len(data['binary_distances'])==sum(r['candidates']for r in record['binary_results']),'binary witness coverage')
for offset,message,claimed in zip(data['ring_offsets'].tolist(),data['ring_messages'].tolist(),data['ring_distances'].tolist()):
 vector(offset,24,8,'offset')
 c=word(message);d=sum(homogeneous((a-b)%8)for a,b in zip(offset,c))
 require(d==claimed,'ring witness mismatch')
 require(0<d<24,'ring offset not rejected by this witness');ring+=1
M=data['binary_masks'].tolist();Q=data['binary_messages'].tolist();L=data['binary_layers'].tolist();D=data['binary_distances'].tolist()
for mask,message,layer,claimed in zip(M,Q,L,D):
 vector(mask,24,16,'mask');require(type(layer)is int and layer in(0,1),'layer domain')
 parity=sum(m.bit_count()%2 for m in mask);require(0<parity<24,'missing disjointness parity condition')
 c=word(message);d=sum((T[a]^m^layer).bit_count()for a,m in zip(c,mask));require(d==claimed and d>0,'binary witness mismatch')
 if d<24:binary+=1
pairs=np.load(OUT/'pair_witnesses.npz');pair_count=0
require(len(pairs['source_indices'])==sum(d>=24 for d in D),'pair survivor coverage')
require(len(set(pairs['source_indices'].tolist()))==len(pairs['source_indices']),'duplicate pair source')
for index,mask,first,second,layer,claimed in zip(pairs['source_indices'].tolist(),pairs['masks'].tolist(),pairs['first_messages'].tolist(),pairs['second_messages'].tolist(),pairs['layers'].tolist(),pairs['distances'].tolist()):
 require(type(index)is int and 0<=index<len(M),'source index domain')
 vector(mask,24,16,'pair mask');require(type(layer)is int and layer in(0,1),'pair layer domain')
 require(M[index]==mask and D[index]>=24,'pair-source mismatch');a=word(first);b=word(second)
 d=sum((T[x]^T[y]^m^layer).bit_count()for x,y,m in zip(a,b,mask))
 require(0<d==claimed<24,'pair rejection mismatch');pair_count+=1
triple=json.loads((OUT/'triple_results.json').read_text());triple_count=0
proposals=json.loads((OUT/'half_cosets.json').read_text())['unverified_triples'];covered=set()
for p in triple['rejection_witnesses']:
 index=p.get('proposal_index');require(type(index)is int and 0<=index<len(proposals) and index not in covered,'proposal index coverage');covered.add(index)
 proposal=proposals[index];vector(p['functional'],12,2,'functional');require(any(p['functional']),'zero functional');vector(p['offset'],24,8,'D offset')
 require(p['functional']==proposal['functional'],'proposal functional mismatch')
 for offset in proposal['offsets']:vector(offset,24,8,'proposal offset')
 require(len(proposal['offsets'])==2,'proposal width')
 allowed=proposal['offsets']+[[ (a-b)%8 for a,b in zip(*proposal['offsets'])]]
 require(p['offset']in allowed,'witness not bound to proposal pair')
 require(sum(a*b for a,b in zip(p['functional'],p['parent_message']))%2==0,'message not in D')
 c=word(p['parent_message']);d=sum(homogeneous((a-b)%8)for a,b in zip(p['offset'],c))
 require(0<d==p['distance']<24,'D-coset rejection mismatch');triple_count+=1
require(binary+pair_count==len(M),'missing binary rejection')
require(covered==set(range(len(proposals)))and triple_count==triple['input_triples'] and triple['survivors']==0,'missing triple rejection')
result=dict(status='PASS',ring_rejection_witnesses=ring,binary_nearest_point_rejections=binary,binary_two_parent_rejections=pair_count,binary_candidates_total=len(M),three_half_coset_rejections=triple_count,unverified_three_half_coset_survivors=triple['survivors'],scope='Actual positive-distance rejection witnesses only. No lower distance bound is inferred from samples.',elapsed_seconds=time.monotonic()-start,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'independent_witness_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
