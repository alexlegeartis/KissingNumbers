"""Supplement the full proof with valid-delta and deliberately invalid replays.

The full baseline proof is replayed by verify.py. These controls exercise the
delta with that pinned baseline, avoiding repeated unchanged baseline replays.
"""
from contextlib import redirect_stdout
from copy import deepcopy
from io import StringIO
from pathlib import Path
import argparse
import json
import subprocess
import sys
import tempfile
import numpy as np
sys.dont_write_bytecode=True
import verify

HERE=Path(__file__).resolve().parent


def controls():
    c=json.loads((HERE/'certificate.json').read_text())
    with redirect_stdout(StringIO()):positive=verify.verify(replay_baseline=False)
    assert (positive['points_36'],positive['points_37'])==(501822,505638)
    cases=[]
    with tempfile.TemporaryDirectory(prefix='ct_exchange_controls_')as td:
        path=Path(td)/'certificate.json'
        def reject(name,mutation,extra=None):
            d=deepcopy(c);mutation(d);path.write_text(json.dumps(d))
            output=StringIO()
            try:
                with redirect_stdout(output):verify.verify(path,replay_baseline=False)
            except (AssertionError,ValueError,KeyError,IndexError)as exc:
                cases.append(dict(case=name,rejected=True,error_type=type(exc).__name__,**(extra or {})))
            else:raise AssertionError('Invalid certificate accepted: '+name)
        reject('duplicate_added_orbit',lambda d:d['add_phase_orbits'].append(d['add_phase_orbits'][0]))
        source=json.loads((verify.PACKAGE/'completion/certificate.json').read_text())
        G=np.array(source['gram_12'],np.int64);M=np.array(source['sublattice_matrix'],np.int64)
        phi=np.eye(12,dtype=np.int64)-M;R0=np.array(source['root_phase_zero_vectors'],np.int64)
        R=np.stack([R0,R0@phi,R0@phi@phi],axis=1);GG=np.kron(np.eye(3,dtype=np.int64),G)
        def expand(v):
            a,b,d,u,w=v
            return np.array([np.r_[R[a,t],R[b,(u+t)%3],R[d,(w+t)%3]]for t in range(3)])
        restored=expand(c['remove_phase_orbits'][-1])
        new=np.concatenate([expand(v)for v in c['add_phase_orbits']])
        maximum=int((restored@GG@new.T).max());assert maximum>6
        def restore(d):
            d['remove_phase_orbits'].pop();d['points_36']+=3;d['points_37']+=3;d['gain_each']+=3
        reject('restore_a_conflicting_removed_orbit',restore,{'exact_bad_product':maximum,'allowed_product':6,'all_point_norms_squared':12})
        def false_count(d):d['points_36']+=1;d['points_37']+=1;d['gain_each']+=1
        reject('inflate_both_counts',false_count)
        reject('alter_pinned_source_hash',lambda d:d['source_sha256'].__setitem__('phase/exact_lattice.py','0'*64))
    disabled=subprocess.run([sys.executable,'-O',str(HERE/'verify.py')],capture_output=True,text=True)
    assert disabled.returncode!=0 and 'Run without -O' in disabled.stderr
    cases.append({'case':'disabled_assertions','rejected':True,'exit_code':disabled.returncode})
    return {'status':'PASS','positive_delta_counts':[501822,505638],
            'positive_delta_passed':True,'baseline_proof_scope':'Pinned; freshly replayed separately by verify.py',
            'negative_controls':len(cases),'controls':cases}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args();report=controls()
    if args.report:args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
