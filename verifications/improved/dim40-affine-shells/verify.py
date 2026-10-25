"""Check complete package integrity and replay the proof in a temporary directory."""
from pathlib import Path
import json, shutil, sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    expected=results()[0]
    with temporary_source(ROOT) as source:
        compiler=shutil.which('g++')
        if compiler is None and Path('C:/msys64/mingw64/bin/g++.exe').exists():
            compiler='C:/msys64/mingw64/bin/g++.exe'
        if compiler is None: raise SystemExit('A C++17 g++ compiler supporting __int128 is required.')
        replay(source,'verify.py','--compiler',compiler)
        receipt=json.loads((source/'replay/verification.json').read_text())
        assert receipt['status']=='PASS' and receipt['points']==expected['total']
    print(json.dumps(dict(status='PASS',results=[expected])))
