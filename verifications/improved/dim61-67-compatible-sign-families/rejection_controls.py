"""Semantic corruptions, including explicit forbidden-inner-product witnesses."""
from copy import deepcopy
from pathlib import Path
import json,subprocess,sys,tempfile
import numpy as np
from verify_family import SOURCE,check,code,generator

HERE=Path(__file__).resolve().parent

def solve(equations):
    pivots={}
    for row,rhs in equations:
        while row:
            p=row.bit_length()-1
            if p not in pivots:pivots[p]=(row,rhs);break
            a,b=pivots[p];row^=a;rhs^=b
        if not row and rhs:return None
    solution=0
    for p,(row,rhs) in sorted(pivots.items()):
        if (row&solution).bit_count()%2!=rhs:solution|=1<<p
    assert all((row&solution).bit_count()%2==rhs for row,rhs in equations)
    return solution

def witness(a,b):
    coordinates=[[i for i in range(r['support'].bit_length()) if r['support']>>i&1] for r in (a,b)]
    equations=[]
    for k,r in enumerate((a,b)):
        equations += [(h<<(16*k),(h&r['shift']).bit_count()%2) for h in r['checks']]
    overlap=sorted(set(coordinates[0])&set(coordinates[1]));assert len(overlap)==10
    equations += [((1<<coordinates[0].index(c))|(1<<(16+coordinates[1].index(c))),0) for c in overlap]
    value=solve(equations);assert value is not None
    signs=[(value>>(16*k))&65535 for k in range(2)]
    dot=sum(4*(-1 if ((signs[0]>>coordinates[0].index(c))^(signs[1]>>coordinates[1].index(c)))&1 else 1) for c in overlap)
    assert dot==40>32
    return dict(inner_product=dot,limit=32,sign_words=signs,supports=[a['support'],b['support']])

def phase_corruption(original,new_pair):
    cert=deepcopy(original);rows=cert['new'];others=cert['new'] if new_pair else cert['old']
    for i,a in enumerate(rows):
        for j,b in enumerate(others):
            if new_pair and j>=i:continue
            overlap=a['support']&b['support']
            if overlap.bit_count()!=10:continue
            aa,bb=code(a,False),code(b,False)
            common=[h for h in aa['dual']&bb['dual'] if h.bit_count()==8 and h&overlap==h]
            assert len(common)==1;h=common[0];c=(h&-h).bit_length()-1
            local=sum(1 for p in range(c) if b['support']>>p&1)
            b['shift']^=1<<local
            return cert,witness(a,b)
    raise AssertionError('No repair pair found')

def scalar_input(cert,path):
    with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False)as z:old=list(map(int,z['words']))
    rows=generator();removed=sum(1<<i for i in cert['removed_coordinates'])
    lines=[f'{len(old)} {len(cert["new"])} {len(cert["old"])} {len(rows)} {int(cert["odd_total_parity"])} {removed}',
           ' '.join(map(str,rows)),' '.join(map(str,old))]
    for group in (cert['old'],cert['new']):lines += [' '.join(map(str,[r['support'],*r['checks'],r['shift']])) for r in group]
    path.write_text('\n'.join(lines)+'\n')

def run():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    global_cert=json.loads((HERE/'data/global64.json').read_text())
    short=json.loads((HERE/'data/short61_62.json').read_text())
    cases=[]
    duplicate=deepcopy(global_cert);duplicate['new'].append(deepcopy(duplicate['new'][0]));cases.append(('duplicate_support',duplicate,None))
    bad=deepcopy(global_cert);bad['new'][0]['checks'][1]^=1;cases.append(('broken_hamming_check',bad,None))
    odd=deepcopy(short);odd['new'][0]['shift']^=1;cases.append(('lost_odd_parity',odd,None))
    for label,new_pair in (('old_new_phase_flip',False),('new_new_phase_flip',True)):
        bad,forbidden=phase_corruption(global_cert,new_pair);cases.append((label,bad,forbidden))
    report=[]
    with tempfile.TemporaryDirectory(prefix='sign-family-controls-') as directory:
        td=Path(directory);binary=td/'check'
        subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'check_family.cpp'),'-o',str(binary)],check=True)
        for label,cert,forbidden in cases:
            path=td/(label+'.json');path.write_text(json.dumps(cert))
            rejected=False
            try:check(path,cpp=False)
            except (AssertionError,RuntimeError,ValueError):rejected=True
            assert rejected,label
            text=td/(label+'.txt');scalar_input(cert,text)
            result=subprocess.run([str(binary),str(text)],capture_output=True,text=True)
            assert result.returncode!=0 and 'certificate rejected' in result.stderr,label
            report.append(dict(case=label,python_rejected=True,cpp_rejected=True,forbidden_pair=forbidden))
    if (HERE/'data/subcode67.json').exists():
        from subcode67 import code as partial_code,verify67
        original=json.loads((HERE/'data/subcode67.json').read_text())
        bad_count=deepcopy(original);bad_count['new_classes'][0]['count']+=1
        partial_cases=[('subcode67_wrong_class_size',bad_count,None)]
        bad=deepcopy(original);found=False
        for a in bad['new_classes']:
            if found:break
            for b in bad['modified_parent_classes']:
                intersection=a['support']&b['support']
                if intersection.bit_count()!=10:continue
                aa,bb=partial_code(a),partial_code(b,True)
                common=[h for h in aa['dual']&bb['dual'] if h.bit_count()==8 and h&intersection==h]
                if not common:continue
                h=common[0];coordinate=(h&-h).bit_length()-1
                index=sum(1 for p in range(coordinate) if b['support']>>p&1)
                b['shift']^=1<<index
                partial_cases.append(('subcode67_wrong_repair_phase',bad,witness(a,b)));found=True;break
        assert found
        with tempfile.TemporaryDirectory(prefix='partial-sign-controls-') as directory:
            for label,bad,forbidden in partial_cases:
                path=Path(directory)/(label+'.json');path.write_text(json.dumps(bad));rejected=False
                try:verify67(path)
                except (ValueError,AssertionError,RuntimeError):rejected=True
                assert rejected,label
                report.append(dict(case=label,python_rejected=True,forbidden_pair=forbidden))
    entrypoints=[('verify_family','check(__import__("pathlib").Path("unused"))'),
                 ('geometry_counts','verify_geometry()'),('section63','verify_section(__import__("pathlib").Path("unused"))'),
                 ('expanded_geometry','check_expanded("unused",66)'),('replay_parents','verify_parents()'),
                 ('baseline67','verify_baseline67()'),('subcode67','verify67("unused")')]
    for name,call in entrypoints:
        result=subprocess.run([sys.executable,'-O','-c',f'from {name} import *; {call}'],cwd=HERE,capture_output=True,text=True)
        assert result.returncode!=0 and 'Run without -O' in result.stderr,name
        report.append(dict(case=name+'_optimized_python',python_rejected=True))
    from verify_subcode67_independent import rejection_controls as independent_controls
    direct=independent_controls(HERE/'data/subcode67.json')
    assert direct['status']=='PASS'
    return dict(status='PASS',controls=len(report)+direct['control_count'],cases=report,
                independent_subcode67_controls=direct)

if __name__=='__main__':print(json.dumps(run(),indent=2))
