"""Exact independent cross-layer audit; baseline replay is a separate prerequisite."""
from pathlib import Path
import json,numpy as np,importlib.util,hashlib
P=Path(__file__).parent;S=P.parent;import argparse
parser=argparse.ArgumentParser();parser.add_argument('--certificate',type=Path,default=P/'certificate.json');parser.add_argument('--points',type=Path,default=P/'new_points.npz');args=parser.parse_args();c=json.loads(args.certificate.read_text());old=json.loads((S/'refill/certificate.json').read_text());assert __debug__;assert c['format']=='orthogonal_hamming_38_v1' and c['dimension']==38 and c['new_point_count']==512 and c['previous_size']==611522;assert hashlib.sha256((S/'refill/certificate.json').read_bytes()).hexdigest()==c['baseline_certificate_sha256']
for name,digest in c['source_sha256'].items():assert hashlib.sha256((S/'source/source'/name).read_bytes()).hexdigest()==digest
assert hashlib.sha256(args.points.read_bytes()).hexdigest()==c['new_points_sha256'];spec=importlib.util.spec_from_file_location('independent_baseline_coordinates',S/'refill/verify_coordinates.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);A,words=m.reconstruct();m.member(A,words);assert A.shape==(98280,24)
roots=np.array(c['body_roots'],np.int64);assert roots.shape==(8,24) and np.array_equal(roots@roots.T,32*np.eye(8,dtype=np.int64));ids=np.array(c['body_line_indices']);assert len(set(ids))==8 and np.array_equal(A[ids],roots);gridids=np.array(old['body_line_indices']);gridowners=set(gridids);removed=set(c['removed_cap_body_lines']);assert len(gridowners)==14 and removed==gridowners|set(ids) and len(removed)==22
z=np.load(args.points);B=z['body_integer'].astype(np.int64);Y=z['head_integer'].astype(np.int64);meta=z['metadata'];assert meta.shape==(512,4) and np.issubdtype(meta.dtype,np.integer);assert np.all((meta[:,0]>=0)&(meta[:,0]<2)) and np.all((meta[:,1]>=0)&(meta[:,1]<128)) and np.all((meta[:,2]>=0)&(meta[:,2]<8)) and np.all(np.isin(meta[:,3],[-1,1]));assert B.shape==(512,24) and Y.shape==(512,14) and len(set(map(tuple,np.c_[B,Y])))==512
expected=set()
for half in range(2):
 for word in range(128):
  syndrome=0
  for k in range(7):
   if word>>k&1:syndrome^=k+1
  y=[0]*14
  for k in range(7):y[7*half+k]=1-2*((word>>k)&1)
  for sign in [-1,1]:expected.add(tuple(list(sign*roots[syndrome])+y))
assert set(map(tuple,np.c_[B,Y]))==expected;assert np.all(np.sum(B*B,axis=1)==32) and np.all(np.sum(Y*Y,axis=1)==7)
# Multiplying all new-new inner products by 21 leaves exact integers.
G=28*(B@B.T)+256*(Y@Y.T);assert np.all(np.diag(G)==2688);np.fill_diagonal(G,-10**9);new_max=int(G.max());assert new_max<=1344
# New versus Fano grid: body products vanish, remaining radical is squared.
assert np.max(np.abs(roots@A[gridids].T))==0;Z=np.array(old['doubled_head_lines'],np.int64);grid_head_max=int(np.abs(Y@Z.T).max());assert 1024*grid_head_max**2<=4096*21
# New versus every retained cap: different minimal Leech lines have |dot|<=16.
retained=np.array([i for i in range(len(A)) if i not in removed]);cap_body_max=int(np.abs(roots@A[retained].T).max());assert cap_body_max<=16;W=np.loadtxt(S/'source/source/dim14_directions.txt',delimiter=',',dtype=np.int64);cap_head_max=int(np.abs(Y@W.T).max());assert cap_head_max<=4
# 2/sqrt(3)+4sqrt(2/21)<=sqrt(6): both squarings preserve signs.
assert 6*7-20==22 and 256*2*49<=22**2*63
# New versus all retained flat axes and coordinate axes.
N=np.array(old['rotation_numerator'],np.int64);axis=W[old['axis_keep']]@N.T;axis_max=int(np.abs(Y@axis.T).max());assert axis_max**2<=21;assert 256<=32*21
# New versus both signs of every equator line, without using a floating bound.
E=np.loadtxt(S/'source/source/equator_lines_6855.txt',dtype=np.int64);eq_max=int(np.abs(roots@E.T).max());assert eq_max**2<=1152
caps=6*(len(A)-len(removed));equator=2*len(E);axes=len(axis)+2*len(old['additional_coordinate_axes']);grid=4*len(gridids);total=caps+equator+axes+grid+len(B);assert total==c['claimed_size']==611638
report={'all_new_and_cross_layer_checks_passed':True,'dimension':38,'lower_bound':total,'gain':total-c['previous_size'],'counts':{'retained_caps':caps,'equator':equator,'axes':axes,'retained_fano_grid':grid,'new_hamming_points':len(B)},'removed_cap_lines':len(removed),'restored_cap_lines':len(set(old['removed_cap_body_lines'])-removed),'newly_removed_cap_lines':len(removed-set(old['removed_cap_body_lines'])),'new_new_pairs':len(B)*(len(B)-1)//2,'maximum_new_new_numerator_over_21':new_max,'new_grid_max_head_integer_product':grid_head_max,'new_cap_max_body_product':cap_body_max,'new_cap_max_head_product':cap_head_max,'new_axis_max_head_product':axis_max,'new_equator_max_body_product':eq_max,'certificate_sha256':hashlib.sha256(args.certificate.read_bytes()).hexdigest(),'arithmetic':'Exact int64 matrix products and integer comparisons after symbolic squaring. Baseline source proof must also be replayed.'};(P/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
