"""Optimize complete E8 head-frame replacements on one 27D tail axis."""
import argparse,itertools,json,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from image_search import above
sys.path.append('/workspace/work/deps')
from ortools.sat.python import cp_model

def canon_frame(F):
    F=F.copy();sg=np.where(F[np.arange(len(F)),np.argmax(F!=0,axis=1)]>0,1,-1);F*=sg[:,None]
    return F[np.lexsort(F.T[::-1])]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--axis',type=int,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=900);ap.add_argument('--seed',type=int,default=641)
    a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=Path('/workspace/work');z=np.load(r/'normal_orbits27_v1/data.npz')
    N,J=z['N'],z['J'];I=np.r_[N,[24,25,26]];H,colors,P,old=z['H'],z['colors'],z['parents'],z['old'];QB,QT,QP=z['QB'],z['QT'],z['quad_parents'];T=z['tail_palette']
    O=2*np.load(r/'PackingStar/25D-31D/new_bounds_20261111/D27_209862_coordinates.npz')['points'];outside=np.any(O[:,J],axis=1);L=np.load(r/'leech.npy');lk={v.tobytes():i for i,v in enumerate(L)}
    original=[]
    for axis in range(3):
        ids=np.flatnonzero(z['old_head']&(z['HT'][:,axis]!=0));assert len(ids)==8;original.append(canon_frame(z['HA'][ids]))
    allheads=np.unique(z['HA'],axis=0);headkeys={v.tobytes() for v in allheads};groups=defaultdict(list)
    for h in allheads:
        rr=h[N]%8;groups[(min(tuple(rr),tuple((rr+4)%8)),int(h[N].sum()%16))].append(h)
    assert len(groups)==135 and all(len(v)==8 for v in groups.values())
    frames=[];framekeys=set()
    def addframe(F):
        F=canon_frame(F).astype(np.int16);key=F.tobytes()
        if key in framekeys:return
        assert np.array_equal(F.astype(int)@F.astype(int).T,64*np.eye(8,dtype=int))
        framekeys.add(key);frames.append(F)
    addframe(original[a.axis])
    for group in groups.values():addframe(np.array(group))
    matrices=[np.ones((4,4),int)-2*np.eye(4,dtype=int),np.array([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]])]
    for ids in itertools.combinations(range(8),4):
        for M in matrices:
            F=original[a.axis].copy();num=M@F[list(ids)].astype(int)
            if np.any(num%2):continue
            new=(num//2).astype(np.int16);new=canon_frame(new)
            if any(h.tobytes() not in headkeys for h in new):continue
            F[list(ids)]=new;addframe(F)
    frames=np.array(frames,np.int16);nf=len(frames);original_keys={v.tobytes() for v in original[a.axis]}
    novelty=np.array([8-sum(h.tobytes() in original_keys for h in frame) for frame in frames],int)
    keep=np.ones(len(H),bool);qkeep=np.ones(len(QB),bool)
    for axis in range(3):
        if axis==a.axis:continue
        hh=np.max(np.abs(H.astype(np.float32)@original[axis].astype(np.float32).T),axis=1)
        qq=np.max(np.abs(QB.astype(np.float32)@original[axis].astype(np.float32).T),axis=1)
        keep&=(T[colors,axis]==0)|(hh<=32);qkeep&=(QT[:,axis]==0)|(qq<=16)
    H,colors,P,old=H[keep],colors[keep],P[keep],old[keep];QB,QT,QP=QB[qkeep],QT[qkeep],QP[qkeep]
    hd=np.max(np.abs(H.astype(np.float32)@frames.reshape(-1,24).astype(np.float32).T).reshape(len(H),nf,8),axis=2)
    qd=np.max(np.abs(QB.astype(np.float32)@frames.reshape(-1,24).astype(np.float32).T).reshape(len(QB),nf,8),axis=2)
    ha=(hd<=32)|(T[colors,a.axis,None]==0);qa=(qd<=16)|(QT[:,a.axis,None]==0)
    keep=ha.any(1);qkeep=qa.any(1);H,colors,P,old,ha=H[keep],colors[keep],P[keep],old[keep],ha[keep];QB,QT,QP,qa=QB[qkeep],QT[qkeep],QP[qkeep],qa[qkeep]
    assert old.sum()==42
    roots=np.flatnonzero(~np.any(L[:,J],axis=1));rootlines=roots[L[roots,np.argmax(L[roots]!=0,axis=1)]>0];assert len(rootlines)==120;where={int(v):i for i,v in enumerate(rootlines)}
    model=cp_model.CpModel();x=[model.new_bool_var(f'e{i}') for i in range(len(H))];q=[model.new_bool_var(f'q{i}') for i in range(len(QB))];f=[model.new_bool_var(f'f{i}') for i in range(nf)]
    used=[model.new_bool_var(f'r{i}') for i in rootlines];double=[model.new_bool_var(f'd{i}') for i in rootlines];model.add_exactly_one(f)
    incident,qowner=[[] for _ in rootlines],[[] for _ in rootlines]
    for i,pp in enumerate(P):
        for p in pp:incident[where[int(p)]].append(i)
    for i,p in enumerate(QP):qowner[where[int(p)]].append(i)
    for j,es in enumerate(incident):
        deg=sum(x[i] for i in es);qq=sum(q[i] for i in qowner[j]);model.add(deg+2*qq<=2*used[j]);model.add(used[j]<=deg+qq);model.add(double[j]==deg+2*qq-used[j])
        for c in range(6):
            same=[x[i] for i in es if colors[i]==c]
            if len(same)>1:model.add_at_most_one(same)
        dd=sum(bool(old[i]) for i in es);model.add_hint(used[j],dd>0);model.add_hint(double[j],max(0,dd-1))
    bodygroups=defaultdict(list);quadgroups=defaultdict(list)
    for i,h in enumerate(H):bodygroups[h.tobytes()].append(i)
    for ids in bodygroups.values():
        if len(ids)>1:model.add_at_most_one([x[i] for i in ids])
        affected=[i for i in ids if T[colors[i],a.axis]!=0]
        if affected:
            bad=np.flatnonzero(~ha[affected[0]])
            if len(bad):model.add_at_most_one([x[i] for i in affected]+[f[int(j)] for j in bad])
    for i,p in enumerate(QP):quadgroups[int(p)].append(i)
    for ids in quadgroups.values():
        affected=[i for i in ids if QT[i,a.axis]!=0]
        if affected:
            bad=np.flatnonzero(~qa[affected[0]])
            if len(bad):model.add_at_most_one([q[i] for i in affected]+[f[int(j)] for j in bad])
    cost=np.abs(T.astype(int)@T.astype(int).T);HF=H.astype(np.float32);QF=QB.astype(np.float32);edges=0
    for at in range(0,len(H),128):
        bad=np.abs(HF[at:at+128]@HF.T)+cost[colors[at:at+128,None],colors[None,:]]>64
        for i,j in zip(*np.nonzero(bad)):
            i+=at
            if i<j:model.add_at_most_one(x[int(i)],x[int(j)]);edges+=1
    bad=above(np.abs(T[colors].astype(np.float32)@QT.astype(np.float32).T),np.abs(HF@QF.T))
    for i,j in zip(*np.nonzero(bad)):model.add_at_most_one(x[int(i)],q[int(j)]);edges+=1
    bad=3*np.abs(QF@QF.T)+np.abs(QT.astype(np.float32)@QT.astype(np.float32).T)>64
    for i,j in zip(*np.nonzero(np.triu(bad,1))):model.add_at_most_one(q[int(i)],q[int(j)]);edges+=1
    D=sum(double);model.add(D<=37).only_enforce_if(f[0]);model.maximize(16*D+sum(int(c)*v for c,v in zip(novelty,f)))
    for i,v in enumerate(x):model.add_hint(v,bool(old[i]))
    for v in q:model.add_hint(v,False)
    for i,v in enumerate(f):model.add_hint(v,i==0)
    def reconstruct(ee,qq,fr):
        removed=np.unique(np.r_[P[ee].ravel(),QP[qq]]);ban=set(map(int,removed));ban.update(lk[(-L[p]).tobytes()] for p in removed)
        rootkeep=np.array([int(p) not in ban for p in roots]);LC=np.c_[2*L[roots[rootkeep]],np.zeros((int(rootkeep.sum()),3),np.int16)]
        mids=np.array([np.r_[s*H[j],t*T[colors[j]]] for j in ee for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,27)
        heads=[]
        for axis in range(3):
            for h in (frames[fr] if axis==a.axis else original[axis]):
                for s,t in itertools.product((-1,1),repeat=2):
                    v=np.r_[s*h,np.zeros(3,np.int16)];v[24+axis]=8*t;heads.append(v)
        QA=np.array([np.r_[np.zeros(24,np.int16),t*QT[j]] for j in qq for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,27)
        QBB=np.array([np.r_[s*QB[j],np.zeros(3,np.int16)] for j in qq for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,27)
        pure=O[~np.any(O[:,:24],axis=1)];AA=np.r_[LC,pure,mids,np.array(heads,np.int16),QA];BB=np.r_[np.zeros_like(LC),np.zeros_like(pure),np.zeros_like(mids),np.zeros_like(np.array(heads,np.int16)),QBB]
        return np.r_[O[outside],AA],np.r_[np.zeros_like(O[outside]),BB]
    print('MODEL',json.dumps(dict(axis=a.axis,frames=nf,bodies=len(H),quads=len(QB),edges=edges,seconds=time.monotonic()-start)),flush=True)
    flags={'completion':False,'positive':False};neutrals=[]
    class Incumbent(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            d=int(self.value(D));fr=next(i for i,v in enumerate(f) if self.boolean_value(v));print('INCUMBENT',d,'normal',348+2*d,'frame',fr,'seconds',self.wall_time,flush=True)
            if d<37 or fr==0:return
            ee=np.flatnonzero([self.boolean_value(v) for v in x]);qq=np.flatnonzero([self.boolean_value(v) for v in q]);RA,RB=reconstruct(ee,qq,fr)
            if d>=38:
                flags['positive']=True;self.stop_search();return
            from augment_neutral27 import augment
            if augment(O,RA,RB,a.out,a.seed+len(neutrals)):flags['completion']=True;self.stop_search()
            neutrals.append(dict(frame=fr,D=d))
    solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.random_seed=a.seed;solver.parameters.linearization_level=2;solver.parameters.log_search_progress=True
    st=solver.solve(model,Incumbent());report=dict(axis=a.axis,status=solver.status_name(st),frames=nf,bodies=len(H),quads=len(QB),neutrals=neutrals,seconds=time.monotonic()-start)
    if flags['completion']:return
    if flags['positive']:
        ee=np.flatnonzero([solver.boolean_value(v) for v in x]);qq=np.flatnonzero([solver.boolean_value(v) for v in q]);fr=next(i for i,v in enumerate(f) if solver.boolean_value(v));RA,RB=reconstruct(ee,qq,fr)
        assert len(RA)>209862 and np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1))
        sourcekeys={v.tobytes() for v in O};newkeys={v.tobytes() for v,b in zip(RA,RB) if not b.any()};removed=np.flatnonzero([v.tobytes() not in newkeys for v in O]);added=np.flatnonzero([b.any() or v.tobytes() not in sourcekeys for v,b in zip(RA,RB)])
        assert len({np.r_[v,b].tobytes() for v,b in zip(RA,RB)})==len(RA)
        A64,B64=RA.astype(np.int64),RB.astype(np.int64)
        for j in added:
            fail=above(A64@A64[j]+3*B64@B64[j],A64@B64[j]+B64@A64[j]);fail[j]=False;assert not fail.any(),('changed pair',j,np.flatnonzero(fail)[:8])
        np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=RA[added],added_B=RB[added],removed_indices=removed)
        report.update(count=len(RA),gain=len(RA)-209862,added=len(added),removed=len(removed));print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
