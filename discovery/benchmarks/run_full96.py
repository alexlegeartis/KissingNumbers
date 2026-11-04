"""Replay both existing full 88--96D layered proofs; not a generic layer search.

This bridge runs the legacy checkers afresh. Comparing saved JSON alone is not
mathematical verification. No record registry or old receipt is modified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/'verifications/improved/dim88-96-z8-golay'
EXPECTED_96=142162362560


def _compare_reports(a,b):
    if a.get('status')!='PASS' or b.get('status')!='PASS':
        raise ValueError('Both freshly executed checkers must pass')
    if a.get('dimensions')!=b.get('dimensions') or set(a.get('dimensions',{}))!={str(d) for d in range(88,97)}:
        raise ValueError('Independent dimension reports disagree or are incomplete')
    row=a['dimensions']['96'];terms=row.get('terms',{})
    if set(terms)!={'dense','transversal24','support24','support6','pure'}:
        raise ValueError('Full96 requires all five layers')
    if any(type(v) is not int or v<=0 for v in terms.values()):
        raise ValueError('Layer cardinalities must be positive integers')
    if sum(terms.values())!=row.get('total') or row['total']!=EXPECTED_96:
        raise ValueError('Full96 control cardinality mismatch')
    if a.get('support_sha256')!=b.get('support_sha256'):
        raise ValueError('Checkers used different support data')
    return row


def run(output_directory,seconds=120):
    if not __debug__:
        raise ValueError('Legacy proofs require assertions; run without -O')
    if type(seconds) is not int or seconds<=0:raise ValueError('Positive per-checker seconds required')
    out=Path(output_directory);out.mkdir(parents=True,exist_ok=False)
    timings={};reports=[];env=dict(os.environ);env['OPENBLAS_NUM_THREADS']='1'
    for name in ('replay','verify_independent'):
        target=out/(name+'.json');start=time.monotonic()
        with (out/(name+'.log')).open('x') as stream:
            subprocess.run([sys.executable,str(PACKAGE/'source'/(name+'.py')),'--output',str(target.resolve())],
                           stdout=stream,stderr=subprocess.STDOUT,env=env,check=True,timeout=seconds)
        timings[name]=time.monotonic()-start
        reports.append(json.loads(target.read_text()))
    row=_compare_reports(*reports)
    paths=[PACKAGE/'source/replay.py',PACKAGE/'source/verify_independent.py',PACKAGE/'data/supports.npz',PACKAGE/'data/z8_generator.json']
    result={'status':'known_full_construction_replayed','dimension':96,'cardinality':row['total'],
            'terms':row['terms'],'per_checker_seconds':timings,'per_checker_time_limit_seconds':seconds,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'scope':'Two existing compressed full-construction checkers, run afresh. Specialized benchmark bridge, not a new universal proof backend or a new record.'}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',required=True);p.add_argument('--per-checker-seconds',type=int,default=120)
    args=p.parse_args();print(json.dumps(run(args.outdir,args.per_checker_seconds),indent=2))
