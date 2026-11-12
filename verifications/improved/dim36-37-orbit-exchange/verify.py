"""Direct integer coordinate check of every affected pair in a phase exchange.

The existing triangle-switch certificate supplies all unchanged pair bounds.
This script imports neither the candidate-admission tables nor search code.
It freshly enumerates the short lattice shells and reconstructs the entire base.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import time
import shutil
import subprocess
import tempfile
import numpy as np

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent/'dim36-37-coxeter-todd'
sys.dont_write_bytecode=True
PINNED_SOURCE={'completion/certificate.json': 'b37691e110c73f1100112aea1a4dc39e2a911559a7a49c67c88bb303f6cd8760', 'completion/verify.py': '4c36c65de1383b5a3db81dc58e51355256fb7ede9afedef269154001bc117db3', 'phase/exact_lattice.py': '5efd39c031287f0f2e3508a4336a715817d7b3786debabd8ea339c95dd5df2d1', 'triangle_switch/certificate.json': '09e5fcfd11cf09a11d68cac48c923928f26591245923921d20a5e7ee8186d16a', 'triangle_switch/proposal.json': '6a88389c8adcdff44b9cd8ef09edaa0bbc9c98af8db4879296c4d4f9ef25356d', 'triangle_switch/switches.json': '54096729335bfcffb792dc06c3d807ad53ceaea939de2b9acbda4f7791729774', 'triangle_switch/verify.py': '48acf68fbafeda8f829f1aea51b606333c20b6fbd03adb93e54c7d024c994ba7'}
sys.path.insert(0,str(PACKAGE/'phase'))


def verify(certificate_path=None,replay_baseline=True):
    if not __debug__:raise RuntimeError('Run without -O')
    certificate_path=Path(certificate_path) if certificate_path else HERE/'certificate.json'
    start=time.monotonic()
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    c=json.loads((PACKAGE/'completion/certificate.json').read_text())
    p=json.loads((PACKAGE/'triangle_switch/proposal.json').read_text())
    meta=json.loads((PACKAGE/'triangle_switch/certificate.json').read_text())
    assert p['source_sha256']==sha(PACKAGE/'completion/certificate.json')
    assert meta['proposal_sha256']==sha(PACKAGE/'triangle_switch/proposal.json')
    old=set(map(tuple,c['phase_orbits']))
    old.difference_update(map(tuple,p['remove_phase_orbits']))
    old.update(map(tuple,p['add_phase_orbits']))
    supplied=json.loads(certificate_path.read_text())
    assert supplied['format']=='ct_current_phase_exchange_v1'
    assert supplied['source_sha256']==PINNED_SOURCE
    for name,digest in PINNED_SOURCE.items():assert sha(PACKAGE/name)==digest,name
    from exact_lattice import enumerate_short
    assert supplied['baseline_proposal_sha256']==meta['proposal_sha256']
    assert supplied['baseline_certificate_sha256']==sha(PACKAGE/'triangle_switch/certificate.json')
    for field in ('remove_phase_orbits','add_phase_orbits'):
        values=supplied[field]
        assert all(len(v)==5 and all(type(x)is int for x in v)for v in values)
        assert len(values)==len(set(map(tuple,values)))
    remove=set(map(tuple,supplied['remove_phase_orbits']))
    add=set(map(tuple,supplied['add_phase_orbits']))
    assert remove<=old and not(add&old)
    rows=np.array(sorted((old-remove)|add),np.int32)
    baseline=None
    if replay_baseline:
        with tempfile.TemporaryDirectory(prefix='ct_orbit_exchange_baseline_')as td:
            temp=Path(td)
            for name in PINNED_SOURCE:
                target=temp/name;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(PACKAGE/name,target)
            subprocess.run([sys.executable,'-B',str(temp/'triangle_switch/verify.py')],check=True,capture_output=True,text=True)
            baseline=json.loads((temp/'triangle_switch/verification.json').read_text())
            assert baseline['all_checks_pass']
            assert (baseline['points_36'],baseline['points_37'])==(501813,505629)
            assert baseline['proposal_sha256']==meta['proposal_sha256']
    assert np.all((rows[:,:3]>=0)&(rows[:,:3]<252))
    assert np.all((rows[:,3:]>=0)&(rows[:,3:]<3))
    candidate=set(map(tuple,rows));assert len(candidate)==len(rows)
    removed=sorted(old-candidate);added=sorted(candidate-old)
    assert len(added)>len(removed)>0
    G=np.array(c['gram_12'],np.int64)
    M=np.array(c['sublattice_matrix'],np.int64)
    GG=np.kron(np.eye(3,dtype=np.int64),G)
    phi=np.eye(12,dtype=np.int64)-M
    R0=np.array(c['root_phase_zero_vectors'],np.int64)
    R=np.stack([R0,R0@phi,R0@phi@phi],axis=1)
    def expand(orbits):
        v=np.array(orbits,np.int64)
        phase=np.c_[np.zeros(len(v),np.int64),v[:,3:]]
        return np.stack([R[v[:,:3],(phase+t)%3]for t in range(3)],axis=1).reshape(-1,36)
    E=expand(sorted(candidate));N=expand(added)
    assert len(np.unique(E,axis=0))==3*len(candidate)
    assert np.all(np.sum((E@GG)*E,axis=1)==12)
    shell=list(enumerate_short(G.tolist()))
    V=np.array([v for v,n in shell],np.int64)
    norms=np.array([n for v,n in shell],np.int64)
    first,second=np.array(c['label_first']),np.array(c['label_second'])
    powers=3**np.arange(6)
    labels=((V[:,first]+V[:,second])%3)@powers
    weights=np.full(729,999,np.int64);weights[0]=0
    np.minimum.at(weights,labels,norms)
    cosets=[np.zeros((1,12),np.int64)if h==0 else V[(labels==h)&(norms==weights[h])]
            for h in range(729)]
    pieces=[];axial=V[norms==4]@M
    for j in range(3):
        A=np.zeros((len(axial),3,12),np.int64);A[:,j]=axial;pieces.append(A)
    assert c['base_outer_cosets'][0]==[0,0,0]
    for row in c['base_outer_cosets'][1:]:
        pieces.append(np.array(list(itertools.product(*(cosets[h]for h in row))),np.int64))
    B=np.concatenate(pieces).reshape(-1,36)
    assert B.shape==(403704,36)
    switches=json.loads((PACKAGE/'triangle_switch/switches.json').read_text())
    assert sorted(s['triangle']for s in switches)==sorted(p['flipped_triangles'])
    Y=np.unique(V[norms==8],axis=0)
    Ylabels=((Y[:,first]+Y[:,second])%3)@powers
    rootlabels=((R0[:,first]+R0[:,second])%3)@powers
    rootmap={int(h):i for i,h in enumerate(rootlabels)}
    owners=[]
    for row in np.array(c['base_outer_cosets'],np.int64):
        if sorted(weights[row])!=[0,4,8]:continue
        i=int(np.flatnonzero(weights[row]==4)[0]);j=int(np.flatnonzero(weights[row]==8)[0]);k=3-i-j
        yi=np.flatnonzero(Ylabels==row[j]);assert len(yi)==81
        owners.append((i,j,k,rootmap[int(row[i])],yi))
    assert len(owners)==756
    O=[];C=[]
    for s in switches:
        a=s['root_index'];i=s['old_root_block'];j=s['norm8_block'];k=s['new_negative_root_block']
        assert sorted([i,j,k])==[0,1,2]
        y=np.array(s['y'],np.int64);assert int(y@G@y)==8
        owner,offset=divmod(s['triangle'],81)
        oi,oj,ok,oa,yi=owners[owner]
        assert (i,j,k,a)==(oi,oj,ok,oa)
        assert (s['owner'],s['Y_offset'])==(owner,offset)
        assert np.array_equal(y,Y[yi[offset]])
        o=np.zeros((3,3,12),np.int64);n=o.copy()
        o[:,i]=R[a];o[:,j]=y;n[:,k]=-R[a];n[:,j]=y
        O.extend(o.reshape(3,36));C.extend(n.reshape(3,36))
    O=np.array(O,np.int64);C=np.array(C,np.int64)
    def keys(a):return np.ascontiguousarray(a).view(np.dtype((np.void,a.shape[1]*a.dtype.itemsize))).ravel()
    killed=np.isin(keys(B),keys(O));assert int(killed.sum())==len(O)==270
    B=np.vstack([B[~killed],C]);assert len(B)==403704
    # Products are int64 sums of 36 integer terms; bound them before arithmetic.
    NG=N@GG
    assert 36*int(np.abs(NG).max())*max(int(np.abs(B).max()),int(np.abs(E).max()))<2**63
    categories={}
    def check(name,left,right,limit,selfkeys=False):
        lm=left@GG;maximum=-999;count=0
        lookup={tuple(v):i for i,v in enumerate(right)}if selfkeys else None
        for a in range(0,len(right),8192):
            dot=lm@right[a:a+8192].T
            if selfkeys:
                for i,v in enumerate(left):
                    j=lookup[tuple(v)]
                    if a<=j<a+len(right[a:a+8192]):dot[i,j-a]=-999
            maximum=max(maximum,int(dot.max()));assert dot.max()<=limit
            count+=dot.size
        if selfkeys:count-=len(left)
        categories[name]={'ordered_products':count,'maximum':maximum,'limit':limit}
        print(name,categories[name],flush=True)
    check('added_vs_current_base',N,B,6)
    check('added_vs_all_final_phase_points',N,E,6,True)
    lift=np.array(c['lift_phase_orbits'],np.int64)
    RX=np.concatenate([R,np.zeros((1,3,12),np.int64)])
    la=np.c_[np.zeros(len(lift),np.int64),lift[:,3:]]
    L=np.stack([RX[lift[:,:3],(la+t)%3]for t in range(3)],axis=1).reshape(-1,36)
    assert len(L)==1908 and np.all(np.sum((L@GG)*L,axis=1)==8)
    check('added_vs_all_latitude_bodies',N,L,6)
    assert len(B)+3*len(old)==meta['points_36']==501813
    points36=len(B)+len(E);points37=points36+2*len(L)
    if supplied is not None:
        assert (points36,points37)==(supplied['points_36'],supplied['points_37'])
    report={'all_affected_pairs_pass':True,'points_36':points36,'points_37':points37,
            'gain_each':points36-meta['points_36'],
            'removed_orbits':len(removed),'added_orbits':len(added),
            'certificate_sha256':sha(certificate_path),
            'short_shell_counts':{str(n):int(np.sum(norms==n))for n in [4,6,8]},
            'categories':categories,'unchanged_pairs':'Inherited from pinned triangle-switch certificate',
            'enumeration':'Fresh rational enumeration of all short lattice shells',
            'arithmetic':'Exact int64 coordinate products with explicit overflow bound',
            'elapsed_seconds':time.monotonic()-start}
    report['baseline_freshly_replayed']=baseline is not None
    report['source_sha256']=PINNED_SOURCE
    if baseline is not None:
        report['baseline_proof']={key:baseline[key]for key in ['points_36','points_37','all_checks_pass','point_pair_checks','source_short_shell_counts']}
    assert supplied['previous_points_36']==501813 and supplied['previous_points_37']==505629
    assert supplied['gain_each']==points36-501813==points37-505629
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,default=HERE/'certificate.json')
    parser.add_argument('--report',type=Path,help='optional JSON receipt path; default replay is read-only')
    args=parser.parse_args()
    report=verify(args.certificate)
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
