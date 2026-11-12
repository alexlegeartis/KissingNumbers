#!/usr/bin/env python3
"""Meaningful corruptions must fail the exact certificate checker."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import shutil
import subprocess
import sys
import verify

HERE=Path(__file__).resolve().parent

def controls():
    if not __debug__: raise RuntimeError('Run without -O.')
    for dim in (68,69): verify.check(HERE/f'certificate{dim}.json')
    cases=[]
    def case(name,dim,kind,mutation):cases.append((name,dim,kind,mutation))
    def field(key,value): return lambda d:d.__setitem__(key,value(d))
    case('changed exact rational mean',69,'certificate',field('exact_value',lambda d:str(F(d['exact_value'])+1)))
    case('changed population ceiling',69,'certificate',field('exact_bound',lambda d:d['exact_bound']+1))
    case('changed final bound',68,'certificate',field('bound',lambda d:d['bound']+1))
    case('changed pole count',69,'certificate',field('with_two_poles',lambda d:d['with_two_poles']+1))
    case('changed dual multiplier',69,'certificate',lambda d:d['dual_exact'].__setitem__(0,str(F(d['dual_exact'][0])+1)))
    case('changed objective coefficient',69,'certificate',lambda d:d['objective'].__setitem__(1,str(F(d['objective'][1])+1)))
    case('removed orbit representative',68,'certificate',lambda d:d['representatives'].pop())
    case('changed moment right side',68,'certificate',lambda d:d['rhs'].__setitem__(2,d['rhs'][2]+1))
    case('invalid head height',69,'certificate',lambda d:d['heights'].__setitem__(0,'3'))
    case('nonzero height in dimension68',68,'certificate',lambda d:d['heights'].__setitem__(0,'1/10'))
    case('changed actual normal vector',69,'normals',lambda d:d['normals'][6].__setitem__(0,d['normals'][6][0]+1))
    case('changed catalogue Gram',69,'embedding',lambda d:d['catalogue_gram'][0].__setitem__(0,d['catalogue_gram'][0][0]+2))
    reports=[]
    for name,dim,kind,mutation in cases:
        with TemporaryDirectory(prefix='gamma72-control-') as tmp:
            directory=Path(tmp)
            for filename in ('embedding.json','d7_embedding.json'):
                shutil.copyfile(HERE/filename,directory/filename)
            c=json.loads((HERE/f'certificate{dim}.json').read_text())
            if kind=='certificate':mutation(c)
            else:
                path=directory/('d7_embedding.json' if kind=='normals' else 'embedding.json')
                d=json.loads(path.read_text());mutation(d);path.write_text(json.dumps(d))
            path=directory/'certificate.json';path.write_text(json.dumps(c))
            rejected=False
            try:verify.check(path,data_dir=directory)
            except (AssertionError,ValueError,KeyError,IndexError,TypeError,StopIteration):rejected=True
            assert rejected,name
            reports.append(dict(control=name,result='REJECTED'))
    process=subprocess.run([sys.executable,'-O',str(HERE/'verify.py')],capture_output=True,text=True)
    assert process.returncode!=0 and 'Run without -O' in process.stderr
    reports.append(dict(control='disabled assertions with -O',result='REJECTED'))
    return dict(status='PASS',baseline_certificates=2,rejection_controls=len(reports),controls=reports)

if __name__=='__main__':
    result=controls()
    print(json.dumps(result,indent=2))
