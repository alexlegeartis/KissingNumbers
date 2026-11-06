"""Both independent routes must reject deliberately corrupted certificates."""
from copy import deepcopy
import importlib.util,json,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    assert __debug__, 'Run without -O'
    verifiers=[]
    for name in ('verify.py','verify_coordinates.py'):
        spec=importlib.util.spec_from_file_location('negative_'+name.replace('.','_'),HERE/name)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        verifiers.append((name,module.verify))
    original=json.loads((HERE/'certificate.json').read_text());cases={}
    bad=deepcopy(original);bad['full_even_supports'][1]=bad['full_even_supports'][0];cases['duplicate support']=bad
    bad=deepcopy(original);bad['heavy_numerators_over_85'][0][0]+=1;cases['wrong heavy norm']=bad
    bad=deepcopy(original);bad['removed_root_lines'].remove([0,2,-1]);cases['restore incompatible root']=bad
    bad=deepcopy(original);bad['points']+=1;cases['false count']=bad
    bad=deepcopy(original);bad['stars'][1]=bad['stars'][0];cases['duplicate star']=bad
    results=[]
    with tempfile.TemporaryDirectory(prefix='rotated34_negative_') as directory:
        path=Path(directory)/'certificate.json'
        for label,bad in cases.items():
            path.write_text(json.dumps(bad))
            for name,verify in verifiers:
                try:verify(path)
                except (AssertionError,ValueError):results.append(dict(case=label,verifier=name,rejected=True))
                else:raise AssertionError((label,name,'accepted corruption'))
    (HERE/'negative_control_results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(dict(all_rejected=True,checks=len(results))))

if __name__=='__main__':main()
