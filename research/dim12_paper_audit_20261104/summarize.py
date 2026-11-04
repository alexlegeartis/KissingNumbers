"""Collect the finite audit experiments without promoting approximate results."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib,json,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def main():
    evidence=ROOT/'evidence';evidence.mkdir(exist_ok=True);report={};best_records=[]
    for name in ['cube49','riesz','bridge_search']:
        rows=json.loads((ROOT/'runs'/name/'summary.json').read_text())
        if name=='cube49':
            b=min(rows,key=lambda r:r['energy']);report[name]=dict(trials=len(rows),best=b,
               certified_successes=sum(r['certified61'] for r in rows))
            assert report[name]['certified_successes']==0
        else:
            assert not any(r.get('exact_certificate',{}).get('passed',False) for r in rows), 'Review a newly certified bound before writing a negative report.'
            b=min(rows,key=lambda r:r['minimax']['max_cos'])
            report[name]=dict(trials=len(rows),best=b)
            groups=['control','placement','rotation','archive'] if name=='riesz' else [False,True]
            group_results={}
            for group in groups:
                data=[r for r in rows if (r['mode'] if name=='riesz' else r['alternative_bridges'])==group]
                best=min(data,key=lambda r:r['minimax']['max_cos'])
                group_results[str(group)]=dict(trials=len(data),best_index=best['index'],best=best['minimax'])
                if name=='bridge_search':best_records.append((name,best))
            report[name]['groups']=group_results
            if name=='riesz':best_records.append((name,b))
    for i,(name,r) in enumerate(best_records):
        filename=f"candidate_{r['index']:04d}.npy" if name=='riesz' else f"bridge_{r['index']:04d}.npy"
        source=ROOT/'runs'/name/filename;target=evidence/f'candidate_{i:02d}.npy';shutil.copyfile(source,target)
        r['saved_file']=str(target.relative_to(ROOT));r['sha256']=hashlib.sha256(target.read_bytes()).hexdigest()
    report.update(verified_lower_bound=841,new_lower_bound_found=False,
                  paper_record_decimal_certificate='decimal841_certificate.json',
                  no_global_impossibility_conclusion=True)
    (ROOT/'experiments.json').write_text(json.dumps(report,indent=2)+'\n')
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    (ROOT/'source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
    print(json.dumps(dict(verified_lower_bound=841,
          cube49_trials=report['cube49']['trials'],full842_trials=report['riesz']['trials']+report['bridge_search']['trials'],
          best_new_max_cos=min(report[x]['best']['minimax']['max_cos'] for x in ['riesz','bridge_search'])),indent=2))


if __name__=='__main__':main()
