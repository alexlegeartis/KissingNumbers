"""Independent SAT encoding of holes and, optionally, unrestricted line labels.

The full encoding permits several true labels at a word: deleting all but one
preserves every constraint, so this is equisatisfiable with one-label geometry.
"""
import argparse, json, time
import numpy as np
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
from model import C, N, DATA, labels_from_holes, verify_labels, verify_component_geometry


def build(a):
    H = (np.load(DATA/'nonlinear_quotient_holes.npy')[a.initial] if a.nonlinear
         else np.load(DATA/'recovery.npz')['holes'][a.initial])
    cnf = CNF(); pool = IDPool(start_from=1025)
    if a.fixed_half is not None:
        from certify_halves import quotient
        _,T=quotient(a.fixed_half)
        for u in np.searchsorted(C,T): cnf.append([int(u)+1])
    if a.half_cuts:
        from certify_halves import verify
        records=json.loads((DATA/'half_certificates.json').read_text())['certificates']
        verify(records)
        kernel=np.load(DATA/'recovery.npz')['holes'][0]
        for r in records:
            half=np.zeros(1024,dtype=bool);half[np.searchsorted(C,r['fixed_subgroup'])]=True
            for fixed in (half,kernel&~half):
                assert a.cap<64
                cnf.append([-int(u)-1 for u in np.flatnonzero(fixed)])
    for u in range(1024):
        cnf.append([u+1]+[int(v)+1 for v in N[u]])
        if a.mode == 'independent':
            for v in N[u]:
                if v > u: cnf.append([-u-1, -int(v)-1])
    x = None
    if a.mode == 'labels':
        x = np.array([[pool.id(('x',u,p)) for p in range(21)] for u in range(1024)])
        for u in range(1024):
            cnf.append([u+1]+list(map(int,x[u])))
            for p in range(21):
                cnf.append([-u-1,-int(x[u,p])])
                cnf.append([-int(x[u,p]), int(N[u,p])+1])
                for v in N[u]:
                    if v > u: cnf.append([-int(x[u,p]),-int(x[v,p])])
    cnf.extend(CardEnc.atmost(list(range(1,1025)), a.cap, vpool=pool,
                             encoding=EncType.cardnetwrk).clauses)
    if a.radius is not None:
        cnf.extend(CardEnc.atmost([-int(u)-1 for u in np.flatnonzero(H)], a.radius,
                                 vpool=pool, encoding=EncType.cardnetwrk).clauses)
    else:
        cnf.append([1])  # translation, not additional hole-set symmetry
    return cnf, H, x


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--name',required=True)
    p.add_argument('--mode',choices=['independent','domination','labels'],default='independent')
    p.add_argument('--cap',type=int,default=63)
    p.add_argument('--radius',type=int)
    p.add_argument('--initial',type=int,default=0)
    p.add_argument('--nonlinear',action='store_true')
    p.add_argument('--hint',action='store_true')
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--solver',default='cadical195')
    p.add_argument('--proof',action='store_true')
    p.add_argument('--fixed-half',type=int)
    p.add_argument('--half-cuts',action='store_true')
    a=p.parse_args();start=time.monotonic()
    cnf,H,x=build(a)
    if a.proof: cnf.to_file(str(DATA/(a.name+'.cnf')))
    print(json.dumps(dict(event='start',run=a.name,variables=cnf.nv,
                          clauses=len(cnf.clauses),seconds=time.monotonic()-start)),flush=True)
    with Solver(name=a.solver,bootstrap_with=cnf,with_proof=a.proof) as s:
        if a.hint:
            phases=[u+1 if H[u] else -u-1 for u in range(1024)]
            if x is not None:
                L=labels_from_holes(H)
                phases += [int(x[u,p])*(1 if L[u]==p else -1)
                           for u in range(1024) for p in range(21)]
            s.set_phases(phases)
        status=None
        while time.monotonic()-start<a.seconds and status is None:
            s.conf_budget(1000)
            status=s.solve_limited()
        report=dict(run=a.name,mode=a.mode,cap=a.cap,radius=a.radius,initial=a.initial,
                    nonlinear=a.nonlinear,hint=a.hint,solver=a.solver,workers=1,
                    fixed_half=a.fixed_half,half_cuts=a.half_cuts,
                    status='SAT' if status else ('UNSAT' if status is False else 'UNKNOWN'),
                    seconds=time.monotonic()-start,variables=cnf.nv,clauses=len(cnf.clauses),
                    stats=s.accum_stats())
        if status:
            positives=set(v for v in s.get_model() if v>0)
            holes=np.array([u+1 in positives for u in range(1024)])
            assert holes.sum()<=a.cap and np.all(holes|holes[N].any(1))
            if a.radius is not None: assert np.sum(H&~holes)<=a.radius
            np.save(DATA/(a.name+'_holes.npy'),holes)
            report['holes']=int(holes.sum())
            if a.mode=='independent': labels=labels_from_holes(holes)
            elif a.mode=='labels':
                labels=np.array([next((p for p in range(21) if x[u,p] in positives),-1)
                                 for u in range(1024)])
                assert np.array_equal(labels<0,holes)
                verify_labels(labels)
            else: labels=None
            if labels is not None:
                np.save(DATA/(a.name+'_labels.npy'),labels)
                report['geometry']=verify_component_geometry(labels)
        if status is False and a.proof:
            proof=s.get_proof()
            (DATA/(a.name+'.drat')).write_text('\n'.join(proof)+'\n')
            report['proof_lines']=len(proof)
            report['proof_checked']=False
        (DATA/(a.name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report),flush=True)
