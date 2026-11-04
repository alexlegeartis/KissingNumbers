"""Preserve compact, reviewable evidence after finite research batches finish."""
from pathlib import Path
import hashlib,json,shutil
from datetime import datetime,timezone
import numpy as np

ROOT=Path(__file__).resolve().parent
BATCHES=['ensemble_fast','continuation_fast','cell_deformations','orbit_search','gram_reflection','orbit_refine','balanced_family','gpu_polish','gpu_final_polish']


def coordinate_path(batch,r,minimax):
    if batch=='ensemble_fast':stem=f"{r['route']}_{r['trial']:04d}"
    elif batch=='continuation_fast':stem=f"{r['method']}_{r['source']}_{r['index']:04d}"
    elif batch=='cell_deformations':stem=f"cells_{r['index']:04d}"
    elif batch=='orbit_search':stem=f"orbit_{r['index']:04d}"
    elif batch=='gram_reflection':return ROOT/'runs'/batch/f"gram_{r['index']:04d}{'' if minimax else '_energy'}.npy"
    elif batch=='orbit_refine':stem=f"pool_{r['index']:04d}"
    elif batch=='balanced_family':stem=f"balanced_{r['index']:04d}"
    elif batch in ['gpu_polish','gpu_final_polish']:stem=f"gpu_{r['index']:04d}"
    return ROOT/'runs'/batch/(stem+('_minimax' if minimax else '')+'.npy')


def main():
    evidence=ROOT/'evidence';evidence.mkdir(exist_ok=True)
    inputs=json.loads((ROOT/'inputs.json').read_text());rows=[];brief=[]
    for batch in BATCHES:
        p=ROOT/'runs'/batch/'summary.json'
        if not p.exists():continue
        a=json.loads(p.read_text())
        for r in a:
            candidate=r.get('minimax',r['best']);use_minimax='minimax' in r
            rows.append(dict(batch=batch,identifier={k:r[k] for k in ['route','method','source','index','trial','seed'] if k in r},
                             energy_best=r['best'],minimax_best=candidate,
                             coordinate_path=str(coordinate_path(batch,r,use_minimax).relative_to(ROOT)),
                             elapsed_seconds=r['elapsed_seconds']))
        energy=min(a,key=lambda r:r['best']['energy'])
        peak=min(a,key=lambda r:r.get('minimax',r['best'])['max_cos'])
        brief.append(dict(batch=batch,completed_trials=len(a),
                          best_energy=energy['best'],best_max_cos=peak.get('minimax',peak['best'])))
        shutil.copyfile(ROOT/'runs'/batch/'config.json',evidence/f'{batch}_config.json')
    native={}
    for batch in ['orbit_tabu','orbit_tabu_weighted','orbit_tabu_partial']:
        p=ROOT/'runs'/batch/'summary.json'
        if p.exists():native[batch]=json.loads(p.read_text())
    # Keep a representative from every completed batch, deduplicating coordinates.
    ranked=sorted(rows,key=lambda r:r['minimax_best']['max_cos'])
    selected=[min([r for r in rows if r['batch']==batch],key=lambda r:r['minimax_best']['max_cos'])
              for batch in BATCHES if any(r['batch']==batch for r in rows)]
    if rows:selected.append(min(rows,key=lambda r:r['energy_best']['energy']))
    saved=[];seen=set()
    for r in selected:
        path=ROOT/r['coordinate_path']
        if not path.exists():continue
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest in seen:continue
        seen.add(digest);target=evidence/f"candidate_{len(saved):02d}.npy";shutil.copyfile(path,target)
        saved.append(dict(file=str(target.relative_to(ROOT)),sha256=digest,**r))
    result=dict(date='2026-11-04',completed_utc=datetime.now(timezone.utc).isoformat(),
                dimension=12,verified_lower_bound=841,new_lower_bound_found=False,
                published842_max_cos=inputs['baseline']['published842']['max_cos'],
                main_construction_trials=sum(r['completed_trials'] for r in brief),batches=brief,
                native_search=native,native_total_iterations=sum(v['iterations'] for v in native.values()),
                saved_candidates=saved,
                exact_whole_cell_result=json.loads((ROOT/'exact_cells_report.json').read_text()))
    hole_report=ROOT/'runs'/'family_holes'/'summary.json'
    if hole_report.exists():
        result['pair_insertion_search']=json.loads(hole_report.read_text())
        shutil.copyfile(hole_report,evidence/'family_holes_summary.json')
    result['gpu_search']={}
    for batch in ['gpu_fire','gpu_fire_final']:
        p=ROOT/'runs'/batch/'summary.json'
        if p.exists():
            result['gpu_search'][batch]=json.loads(p.read_text())
            shutil.copyfile(p,evidence/f'{batch}_summary.json')
    if rows:result['best_final_max_cos']=ranked[0]['minimax_best']['max_cos']
    if any(v['found_97_compatible_orbits'] for v in native.values()):
        raise RuntimeError('A 97-orbit construction needs its separate exact algebraic certificate.')
    for r in ranked:
        if r['minimax_best']['max_cos']>.5:break
        from verify import verify
        certificate=verify(ROOT/r['coordinate_path'])
        if certificate['passed']:
            result['new_lower_bound_found']=True
            result['verified_lower_bound']=max(result['verified_lower_bound'],certificate['n'])
    (ROOT/'campaign.json').write_text(json.dumps(result,indent=2)+'\n')
    (evidence/'trial_metrics.json').write_text(json.dumps(rows,indent=2)+'\n')
    for name in ['orbit_pool/report.json','orbit_pool/weights_report.json','quaternion_cells/report.json']:
        p=ROOT/'runs'/name
        if p.exists():shutil.copyfile(p,evidence/(name.replace('/','_')))
    (ROOT/'source_manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(ROOT.iterdir()) if p.suffix in ['.py','.cpp','.txt'] and p.is_file()},indent=2)+'\n')
    # The declarations above are explicitly fail-closed: a new bound requires an
    # independent exact verification and an intentional status update.
    print(json.dumps({k:result[k] for k in ['main_construction_trials','native_total_iterations','verified_lower_bound','new_lower_bound_found','best_final_max_cos']},indent=2))


if __name__=='__main__':main()
