"""Hash integrity, derived counts, then isolated independent exact proofs."""
from pathlib import Path
import json,sys,os,subprocess
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source
from final import results
SCRIPTS={26:'verify26_201108.py',27:'verify209532.py',28:'verify28_216680.py',29:'verify29_226520_bounded.py',30:'verify30_238796_full_final.py',31:'verify31_246596_full_final.py'}
if __name__=='__main__':
    check_manifest(ROOT)
    expected=results()
    dims=list(map(int,sys.argv[1:])) or list(SCRIPTS)
    assert all(n in SCRIPTS for n in dims)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONHASHSEED='0')
    with temporary_source(ROOT) as source:
        p=source/'27d_independent_verification_20261026'
        for n in dims:
            subprocess.run([sys.executable,'-B',str(p/SCRIPTS[n])],cwd=p,env=env,check=True,timeout=1800)
    print(json.dumps(dict(status='PASS',dimensions=dims,results=[r for r in expected if r['dimension'] in dims])))
