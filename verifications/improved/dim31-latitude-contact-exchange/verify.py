"""Fresh complete geometry, independent direct-integer delta, and file manifest."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile
import numpy as np
from construction import ROOT,delta,reconstruct,require
from verify_geometry import verify,canonical_hash

def manifest():
    path=ROOT/'data/SHA256SUMS.txt'
    expected={name:sha for sha,name in (line.split('  ',1) for line in path.read_text().splitlines() if line)}
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p!=path and p.suffix!='.pyc'}
    require(set(expected)==actual,'Complete manifest coverage')
    for name,sha in expected.items():require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,'File hash: '+name)
    return len(expected)

def changed_pairs(final):
    _,added,_=delta();require(len({v.tobytes() for v in final})==len(final),'Distinct complete coordinates')
    require(np.all(np.sum(final.astype(np.int64)**2,axis=1)==128),'Integer norms')
    lookup={v.tobytes():i for i,v in enumerate(final)};mx=-128
    for row in added:
        dots=np.zeros(len(final),np.int64)
        for j in np.flatnonzero(row):dots+=int(row[j])*final[:,j].astype(np.int64)
        dots[lookup[row.tobytes()]]=-128;require(int(dots.max())<=64,'Independent int64 inserted angle');mx=max(mx,int(dots.max()))
    return {'status':'PASS','dimension':31,'points':len(final),'removed':4,'added':8,'gain':4,'direct_integer_pairs':8*(len(final)-1),'maximum_inserted_product':mx,'canonical_point_sha256':canonical_hash(final)}

def main():
    p=argparse.ArgumentParser();p.add_argument('--write-reports',action='store_true');a=p.parse_args()
    if not a.write_reports:print(json.dumps({'manifest_files':manifest()}),flush=True)
    final=reconstruct();actual=np.load(ROOT/'data/coordinates31_251200.npz',allow_pickle=False)['points'];require(np.array_equal(final,actual),'Reconstruction equals complete witness')
    geometry=verify(ROOT/'data/coordinates31_251200.npz');delta_report=changed_pairs(final)
    for name,report in [('geometry.json',geometry),('changed_pairs.json',delta_report)]:
        path=ROOT/'data'/name
        if a.write_reports:path.write_text(json.dumps(report,indent=2)+'\n')
        else:
            expected=json.loads(path.read_text());require({k:v for k,v in report.items() if k!='seconds'}=={k:v for k,v in expected.items() if k!='seconds'},'Fresh '+name+' agrees')
    print(json.dumps({'status':'PASS','geometry':geometry,'changed_pairs':delta_report}),flush=True)
    # Independently authored checker, explicit portable inputs, fresh full scan.
    try:
        from numpy._core._multiarray_umath import __cpu_features__ as features
    except ImportError:
        try:
            from numpy.core._multiarray_umath import __cpu_features__ as features
        except ImportError:
            features={}
    flags=[] if features.get('AVX512F') and features.get('AVX512BW') else ['--scalar']
    with tempfile.TemporaryDirectory(prefix='latitude31_independent_') as temp:
        command=[sys.executable,str(ROOT/'independent_audit31.py'),
                 str(ROOT/'data/baseline31_251196.npz'),str(ROOT/'data/coordinates31_251200.npz'),
                 '--outdir',temp]+flags
        subprocess.run(command,check=True,cwd=ROOT)
        for name in ('independent_audit31_report.json','independent_full_pairs.json'):
            fresh=json.loads((Path(temp)/name).read_text());path=ROOT/'data'/name
            if a.write_reports:path.write_text(json.dumps(fresh,indent=2)+'\n')
            else:
                old=json.loads(path.read_text());ignored={'seconds','full_integer_scan_seconds','audit_other_seconds','integer_kernel'}
                require({k:v for k,v in fresh.items() if k not in ignored}=={k:v for k,v in old.items() if k not in ignored},'Fresh independent '+name+' agrees')
if __name__=='__main__':main()
