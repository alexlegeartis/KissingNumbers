from pathlib import Path
import numpy as np,json,hashlib
import sys
P=Path(__file__).resolve().parent;W=P.parent;E=P/'medium_families';sys.path.insert(0,str(W));import verify_families as vf
g,residues,mu,tiny=vf.geometry();f=np.load(W/'configuration_43.npz');old=np.load(W/'explicit_points_43.npz');raw=np.load(W/'tiny_bodies.npz')['H'].astype(np.int64);body=np.load(W/'body_complement.npz');cg={'body_gram':body['q'],'body_residues':residues};C=f['classes'];R=f['r3'];Z=f['head_numerators'].astype(object);d=int(f['head_denominator']);OZ=old['head_numerators'].astype(object);od=int(old['head_denominator']);I=old['body_indices'];pure=f['pure_directions'].astype(object);D=json.loads((P/'candidate.json').read_text());den=D['head_denominator'];points=D['points']+D['new_tiny'];FN=R.astype(object)*d*d+np.sum(Z*Z,axis=1);ON=10*od*od+np.sum(OZ*OZ,axis=1);PN=np.sum(pure*pure,axis=1);new=[];checks=0
b=body['b'];q=cg['body_gram'];J=np.array(vf.inverse(q,3),dtype=np.int64);assert np.array_equal(q@J,3*np.eye(40,dtype=np.int64));lookup={tuple(x%3):i for i,x in enumerate(cg['body_residues'])}
medium_body_ids=[]
for k in np.where(R==12)[0]:
 data=np.load(E/f'family_{int(C[k])}.npz');H=data['H'].astype(np.int64);assert H.shape==(495,48) and len(np.unique(H,axis=0))==495;assert np.all(np.sum(H*H,axis=1)==3456);num=H@b.T@J;assert not np.any(num%288);coeff=num//288;assert np.array_equal(coeff@b,H);assert np.all(coeff%3==residues[C[k]]);assert np.array_equal(coeff,data['coefficients']);assert int(data['coset'])==C[k];medium_body_ids.extend(k*495+i for i in range(495))
assert len(set(D['delete']))==len(D['delete']) and set(D['delete'])<=set(medium_body_ids)
assert len(set(D['delete_tiny']))==len(D['delete_tiny']) and all(0<=i<len(I) for i in D['delete_tiny'])
assert den>0
for p in points:
 if 'family' in p:
  j=p['family'];row=p['row'];H=np.load(E/f'family_{int(C[j])}.npz')['H'].astype(np.int64);h=H[row];r=12
 else:
  h=raw[p['body_index']];r=10
 assert int(h@h)==r*288;num=h@b.T@J;assert not np.any(num%288);a=num//288;assert np.array_equal(a@b,h);cl=lookup[tuple(a%3)]
 if 'family' in p:assert cl==C[j]
 z=np.array(p['head_numerators'],dtype=object);n=r*den*den+sum(z*z);assert n==p['normal'];dc=((g['digits'][C]-g['digits'][cl])%3)@g['powers'];beta=r+R-18+g['mu'][dc]//4
 for k in np.where(R!=12)[0]:
  v=int(beta[k])*den*d+2*sum(z*Z[k]);assert v<=0 or v*v<=n*FN[k];checks+=1
 for k in np.where(R==12)[0]:
  other=np.load(E/f'family_{int(C[k])}.npz')['H'].astype(np.int64);dots=other@h;assert not np.any(dots%144);v=(dots//144).astype(object)*den*d+2*sum(z*Z[k]);keep=np.array([k*495+i not in D['delete'] for i in range(len(other))]);assert np.all((v[keep]<=0)|(v[keep]*v[keep]<=n*FN[k]));checks+=int(sum(keep))
 keep_t=np.array([i not in D['delete_tiny'] for i in range(len(I))]);dots=raw[I]@h;assert not np.any(dots%144);v=(dots//144).astype(object)*den*od+2*OZ@z;assert np.all((v[keep_t]<=0)|(v[keep_t]*v[keep_t]<=n*ON[keep_t]));checks+=int(sum(keep_t))
 v=2*pure@z;assert np.all((v<=0)|(v*v<=n*PN));checks+=len(v);new.append((h,z,n))
for i,(h,z,n) in enumerate(new):
 for hh,zz,nn in new[:i]:
  dot=int(h@hh);assert dot%144==0;v=dot//144*den*den+2*sum(z*zz);assert v<=0 or v*v<=n*nn;checks+=1
assert len(new)-len(set(D['delete']))-len(set(D['delete_tiny']))==D['gain']
report={'status':'PASS','old_count':2919799,'new_count':2919799+D['gain'],'new_points':len(new),'deleted_medium_points':len(D['delete']),'deleted_tiny_points':len(D['delete_tiny']),'exact_angle_checks':checks,'candidate_sha256':hashlib.sha256((P/'candidate.json').read_bytes()).hexdigest(),'scope':'Fresh exact checks of new points against every retained explicit medium and tiny point, all pure heads, and universal bounds for retained major families. Published old configuration and population premises remain inherited.'};(P/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
