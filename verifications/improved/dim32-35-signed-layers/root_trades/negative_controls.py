"""Both proof representations must reject corrupted root-star certificates."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import tempfile

HERE=Path(__file__).resolve().parent


def main():
    assert __debug__, 'Run without -O'
    verifiers=[]
    for name in ('verify.py','verify_coordinates.py'):
        spec=importlib.util.spec_from_file_location('negative_'+name.replace('.','_'),HERE/name)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        verifiers.append((name,module.verify))
    checked=0
    with tempfile.TemporaryDirectory(prefix='root_star_negative_') as temporary:
        path=Path(temporary)/'bad.json'
        for dimension in (33,34,35):
            original=json.loads((HERE/('certificate%d.json'%dimension)).read_text())
            cases={}
            bad=deepcopy(original);bad['stars'].append(bad['stars'][0]);bad['points']+=16
            cases['duplicate star']=bad
            bad=deepcopy(original);bad['removed_root_lines'].pop();bad['points']+=2
            cases['undeleted incompatible root']=bad
            bad=deepcopy(original);bad['points']+=1
            cases['false count']=bad
            for label,bad in cases.items():
                path.write_text(json.dumps(bad))
                for name,verify in verifiers:
                    try:
                        verify(path)
                    except (AssertionError,ValueError):
                        checked+=1
                        print('REJECTED:',dimension,label,'by',name,flush=True)
                    else:
                        raise AssertionError((dimension,label,name,'incorrectly accepted'))
    print('PASS:',checked,'independent rejection checks')


if __name__=='__main__':
    main()
