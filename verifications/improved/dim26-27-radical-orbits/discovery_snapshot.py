"""Complete orbit exchanges with variable head frames and radical directions."""
import argparse,itertools,json,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from image_search import above
sys.path.append('/workspace/work/deps')
from ortools.sat.python import cp_model

def positive(v):
    return v if v[np.flatnonzero(v)[0]]>0 else -v

def read_pool(path):
    p=np.load(path)
    A=p['candidate_A'] if 'candidate_A' in p else p['candidates']
    B=p['candidate_B'] if 'candidate_B' in p else np.zeros_like(A)
    den=p['denominators'] if 'denominators' in p else np.full(len(A),int(p['denominator']) if 'denominator' in p else 1)
    cost=np.diff(p['blocker_indptr']) if 'blocker_indptr' in p else np.zeros(len(A),int)
    return p,A,B,den,cost

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('dim',type=int);ap.add_argument('pools',nargs='+',type=Path)
    ap.add_argument('--heads',nargs='*',type=Path,default=[]);ap.add_argument('--quads',nargs='*',type=Path,default=[])
    ap.add_argument('--head-cost',type=int,default=4);ap.add_argument('--quad-cost',type=int,default=6)
    ap.add_argument('--body-cost',type=int,default=2147483647)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=1200)
    ap.add_argument('--max-old-deletions',type=int,default=32);ap.add_argument('--head-deletions',type=int,default=8)
    ap.add_argument('--seed',type=int,default=457);ap.add_argument('--baseline-only',action='store_true')
    a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();root=Path('/workspace/work');d=a.dim
    baseline=201180 if d==26 else 209862
    z=np.load(root/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{baseline}_coordinates.npz')
    O=2*np.r_[z['rational_points'],z['quadratic_A']] if d==26 else 2*z['points']
    OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']] if d==26 else np.zeros_like(O)
    L=np.load(root/'leech.npy');lk={v.tobytes():i for i,v in enumerate(L)}
    anti=np.array([lk[(-v).tobytes()] for v in L])
    pos=L[np.arange(len(L)),np.argmax(L!=0,axis=1)]>0;canon=np.where(pos,np.arange(len(L)),anti)
    g=np.load(root/f'parentgraph{d}/graph.npz');sf=g['source_for_leech'];palette=[]
    for i,j in itertools.combinations(range(d-24),2):
        for s in (-1,1):
            t=np.zeros(d-24,np.int16);t[[i,j]]=4,4*s;palette.append(t)
    T=np.array(palette);tl={t.tobytes():i for i,t in enumerate(T)}
    bodies,parents,original_color={},{},{}
    for h,t,pp in zip(g['bodies'],g['tails'],g['parents']):
        if h[np.flatnonzero(h)[0]]<0:continue
        key=h.tobytes();c=tl[positive(t).tobytes()]
        bodies[key]=(h,{c,c^1});parents[key]=np.sort(canon[pp]);original_color[key]=c
    for path in a.pools:
        p,AA,BB,den,cost=read_pool(path);pp=p['parent_indices'] if 'parent_indices' in p else None
        for i in np.flatnonzero((~np.any(BB,axis=1))&(cost<=a.body_cost)):
            if np.any(AA[i]%int(den[i])):continue
            v=(AA[i]//int(den[i])).astype(np.int16)
            if np.sum(v[:24].astype(int)**2)!=96 or not np.any(v[24:]):continue
            tt=positive(v[24:])
            if tt.tobytes() not in tl:continue
            h=positive(v[:24]);key=h.tobytes();c=tl[tt.tobytes()]
            if key not in bodies:bodies[key]=(h,set())
            bodies[key][1].update((c,c^1))
            if pp is not None:parents[key]=np.sort(canon[pp[i]])
    unknown=[k for k in bodies if k not in parents];LF=L.astype(np.float32)
    for at in range(0,len(unknown),64):
        keys=unknown[at:at+64];dots=np.array([bodies[k][0] for k in keys],np.float32)@LF.T
        for i,key in enumerate(keys):
            pp=np.flatnonzero(dots[i]==48)
            if len(pp)==2 and np.array_equal(L[pp].sum(0),bodies[key][0]):parents[key]=np.sort(canon[pp])
    vals,options=[],[]
    for key,(h,colors) in bodies.items():
        if key not in parents:continue
        group=[]
        for c in sorted(colors):
            group.append(len(vals));vals.append((h,c,parents[key],original_color.get(key,-1)==c))
        options.append(group)
    H=np.array([v[0] for v in vals],np.int16);colors=np.array([v[1] for v in vals],int)
    P=np.array([v[2] for v in vals],int);old=np.array([v[3] for v in vals],bool);origids=np.flatnonzero(old)
    core=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1)
    headgroups,quadgroups={},{}
    for i in np.flatnonzero(core==64):
        h,tt=positive(O[i,:24]),positive(O[i,24:]);headgroups[np.r_[h,tt].tobytes()]=(h,tt,True)
    for i in np.flatnonzero(np.any(OB,axis=1)):
        b,tt=positive(OB[i,:24]),positive(O[i,24:]);quadgroups[np.r_[b,tt].tobytes()]=(b,tt,True)
    for path in a.heads:
        _,AA,BB,den,cost=read_pool(path)
        for i in np.flatnonzero((cost<=a.head_cost)&(~np.any(BB,axis=1))):
            if np.any(AA[i]%int(den[i])):continue
            v=(AA[i]//int(den[i])).astype(np.int16)
            if np.sum(v[:24].astype(int)**2)!=64 or np.sum(v[24:].astype(int)**2)!=64:continue
            h,tt=positive(v[:24]),positive(v[24:]);key=np.r_[h,tt].tobytes()
            headgroups.setdefault(key,(h,tt,False))
    for path in a.quads:
        p,AA,BB,den,cost=read_pool(path)
        if int(p['radicand'])!=3:continue
        for i in np.flatnonzero(cost<=a.quad_cost):
            if den[i]!=1 or np.any(AA[i,:24]) or np.any(BB[i,24:]):continue
            b,tt=positive(BB[i,:24]),positive(AA[i,24:])
            assert b.tobytes() in lk and tt.tobytes() in tl
            quadgroups.setdefault(np.r_[b,tt].tobytes(),(b,tt,False))
    hv,qv=list(headgroups.values()),list(quadgroups.values())
    HA=np.array([v[0] for v in hv],np.int16).reshape(-1,24);HT=np.array([v[1] for v in hv],np.int16).reshape(-1,d-24)
    ho=np.array([v[2] for v in hv],bool)
    QB=np.array([v[0] for v in qv],np.int16).reshape(-1,24);QT=np.array([v[1] for v in qv],np.int16).reshape(-1,d-24)
    qo=np.array([v[2] for v in qv],bool);QP=np.array([canon[lk[v.tobytes()]] for v in QB],int)
    pure=np.flatnonzero(core==0);vertices=np.unique(np.r_[P.ravel(),QP]);where={int(p):i for i,p in enumerate(vertices)}
    incident,qowner=[[] for _ in vertices],[[] for _ in vertices]
    for i,pp in enumerate(P):
        for p in pp:incident[where[int(p)]].append(i)
    for i,p in enumerate(QP):qowner[where[int(p)]].append(i)
    origdeg=np.bincount([where[int(p)] for p in P[old].ravel()],minlength=len(vertices))
    oldquad=np.bincount([where[int(p)] for p in QP[qo]],minlength=len(vertices))
    score=(baseline-len(L)-len(pure))//2
    assert np.sum(origdeg+2*oldquad-(sf[vertices]<0))+2*ho.sum()==score
    def reconstruct(chosen,chosen_heads,chosen_quads):
        removed=np.unique(np.r_[P[chosen].ravel(),QP[chosen_quads]])
        keep=np.ones(len(L),bool);keep[removed]=False;keep[anti[removed]]=False
        LC=np.c_[2*L[keep],np.zeros((int(keep.sum()),d-24),np.int16)]
        C=np.array([np.r_[s*H[j],t*T[colors[j]]] for j in chosen for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,d)
        HC=np.array([np.r_[s*HA[j],t*HT[j]] for j in chosen_heads for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,d)
        QA=np.array([np.r_[np.zeros(24,np.int16),t*QT[j]] for j in chosen_quads for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,d)
        QBB=np.array([np.r_[s*QB[j],np.zeros(d-24,np.int16)] for j in chosen_quads for s,t in itertools.product((-1,1),repeat=2)],np.int16).reshape(-1,d)
        return np.r_[LC,O[pure],C,HC,QA],np.r_[np.zeros_like(LC),OB[pure],np.zeros_like(C),np.zeros_like(HC),QBB]
    RA,RB=reconstruct(origids,np.flatnonzero(ho),np.flatnonzero(qo))
    sourcekeys={np.r_[v,b].tobytes() for v,b in zip(O,OB)}
    assert len(RA)==len(O) and {np.r_[v,b].tobytes() for v,b in zip(RA,RB)}==sourcekeys
    print('EXACT_BASELINE_RECONSTRUCTION',d,len(RA),flush=True)
    if a.baseline_only:return
    for at in range(0,len(HA),64):assert np.max(2*HA[at:at+64].astype(np.float32)@LF.T)<=64
    model=cp_model.CpModel();x=[model.new_bool_var(f'e{i}') for i in range(len(H))]
    heads=[model.new_bool_var(f'h{i}') for i in range(len(HA))];quads=[model.new_bool_var(f'q{i}') for i in range(len(QB))]
    used=[model.new_bool_var(f'r{i}') for i in vertices];double=[model.new_bool_var(f'd{i}') for i in vertices]
    for j,es in enumerate(incident):
        deg=sum(x[i] for i in es);qq=sum(quads[i] for i in qowner[j])
        model.add(deg+2*qq<=2*used[j]);model.add(used[j]<=deg+qq);model.add(double[j]==deg+2*qq-used[j])
        model.add_hint(used[j],bool(sf[vertices[j]]<0));model.add_hint(double[j],int(origdeg[j]+2*oldquad[j]-(sf[vertices[j]]<0)))
        for c in range(len(T)):
            group=[x[i] for i in es if colors[i]==c]
            if len(group)>1:model.add_at_most_one(group)
    for group in options:
        if len(group)>1:model.add_at_most_one([x[i] for i in group])
    for arr,hints in ((x,old),(heads,ho),(quads,qo)):
        for v,b in zip(arr,hints):model.add_hint(v,bool(b))
    if a.max_old_deletions:model.add(sum(1-x[int(i)] for i in origids)<=a.max_old_deletions)
    if a.head_deletions>=0:model.add(sum(1-heads[int(i)] for i in np.flatnonzero(ho))<=a.head_deletions)
    assert np.all(np.count_nonzero(HT,axis=1)==1)
    head_axis=np.argmax(HT!=0,axis=1)
    for axis in range(d-24):
        model.add(sum(heads[int(i)] for i in np.flatnonzero(head_axis==axis))<=24)
    target=sum(double)+2*sum(heads);model.add(target>=score+1)
    HF,HAF,QBF=H.astype(np.float32),HA.astype(np.float32),QB.astype(np.float32);cuts=defaultdict(int)
    def conflicts(leftvars,rightvars,key,batch):
        for at in range(0,len(leftvars),128):
            ii,jj=np.nonzero(batch(at,min(len(leftvars),at+128)))
            for i,j in zip(ii+at,jj):
                if leftvars is rightvars and i>=j:continue
                model.add_at_most_one(leftvars[int(i)],rightvars[int(j)]);cuts[key]+=1
    cost=np.abs(T.astype(np.int64)@T.T)
    conflicts(x,x,'body_body',lambda at,end:np.abs(HF[at:end]@HF.T)+cost[colors[at:end,None],colors[None,:]]>64)
    print('BODY_CONFLICTS',dict(cuts),'seconds',time.monotonic()-start,flush=True)
    conflicts(x,heads,'body_head',lambda at,end:np.abs(HF[at:end]@HAF.T)+np.abs(T[colors[at:end]].astype(np.float32)@HT.astype(np.float32).T)>64)
    conflicts(heads,heads,'head_head',lambda at,end:np.abs(HAF[at:end]@HAF.T)+np.abs(HT[at:end].astype(np.float32)@HT.astype(np.float32).T)>64)
    conflicts(x,quads,'body_quad',lambda at,end:above(np.abs(T[colors[at:end]].astype(np.float32)@QT.astype(np.float32).T),np.abs(HF[at:end]@QBF.T)))
    conflicts(heads,quads,'head_quad',lambda at,end:above(np.abs(HT[at:end].astype(np.float32)@QT.astype(np.float32).T),np.abs(HAF[at:end]@QBF.T)))
    conflicts(quads,quads,'quad_quad',lambda at,end:3*np.abs(QBF[at:end]@QBF.T)+np.abs(QT[at:end].astype(np.float32)@QT.astype(np.float32).T)>64)
    print('MODEL',json.dumps(dict(dimension=d,bodies=len(H),heads=len(HA),quads=len(QB),old_heads=int(ho.sum()),old_quads=int(qo.sum()),vertices=len(vertices),cuts=dict(cuts),target_score=score+1,seconds=time.monotonic()-start)),flush=True)
    class Improve(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            print('FEASIBLE_COMPLETE_MODEL',int(self.value(target)),'seconds',self.wall_time,flush=True);self.stop_search()
    solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=a.seconds
    solver.parameters.random_seed=a.seed;solver.parameters.symmetry_level=0;solver.parameters.linearization_level=2
    solver.parameters.repair_hint=True;solver.parameters.hint_conflict_limit=10000;solver.parameters.log_search_progress=True
    st=solver.solve(model,Improve());report=dict(dimension=d,status=solver.status_name(st),bodies=len(H),heads=len(HA),quads=len(QB),cuts=dict(cuts),seconds=time.monotonic()-start)
    if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        ee=np.flatnonzero([solver.boolean_value(v) for v in x]);hh=np.flatnonzero([solver.boolean_value(v) for v in heads]);qq=np.flatnonzero([solver.boolean_value(v) for v in quads])
        RA,RB=reconstruct(ee,hh,qq)
        assert len(RA)>baseline and len(RA)==len(L)+len(pure)+2*solver.value(target)
        assert np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128)
        assert not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1))
        keys=[np.r_[v,b].tobytes() for v,b in zip(RA,RB)];assert len(set(keys))==len(RA)
        added=np.flatnonzero([k not in sourcekeys for k in keys]);newkeys=set(keys)
        removed=np.flatnonzero([np.r_[v,b].tobytes() not in newkeys for v,b in zip(O,OB)])
        A64,B64=RA.astype(np.int64),RB.astype(np.int64)
        for j in added:
            fail=above(A64@A64[j]+3*B64@B64[j],A64@B64[j]+B64@A64[j]);fail[j]=False
            assert not fail.any(),('changed geometry',j,np.flatnonzero(fail)[:8])
        np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
        np.savez_compressed(a.out/'exchange.npz',added_A=RA[added],added_B=RB[added],removed_indices=removed)
        report.update(count=len(RA),gain=len(RA)-baseline,added=len(added),removed=len(removed));print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)

if __name__=='__main__':main()
