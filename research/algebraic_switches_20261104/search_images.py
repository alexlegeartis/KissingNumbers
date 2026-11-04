"""Whole mixed core/tail isometries of current 27--31D configurations.

Unlike the earlier sampled 31D three-core/one-tail experiment, every source
point is transformed here. A blocker cutoff restricts the search explicitly.
The source configuration is a pinned registered certificate, not assumed to
be a newly discovered construction. Positive outputs require full replay.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'common'))
from golay import golay24
from leech import build

def load_module(path):
    spec=importlib.util.spec_from_file_location('_switch_exchange',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

def maps(n, trials, seed, extended=False):
    rng=np.random.default_rng(seed);code=golay24().astype(np.int64)
    octad=np.array([0,2,4,5,6,10,11,23])
    for i in range(trials):
        M=np.eye(n,dtype=np.int64);denominator=1
        if extended and i%3==0:
            # The cyclic QR automorphisms move the normal octad itself.
            a=int(rng.choice(sorted({j*j%23 for j in range(1,23)})))
            b=int(rng.integers(23));p=[(a*j+b)%23 for j in range(23)]+[23]
            words={tuple(row) for row in code}
            if any(tuple(row[p]) not in words for row in code):raise ValueError('Golay QR automorphism')
            M[:24,:24]=np.eye(24,dtype=np.int64)[p]
            M[:24,:24]*=(1-2*code[int(rng.integers(1,len(code)))])[None,:]
            M[24:,24:]=np.eye(n-24,dtype=np.int64)[rng.permutation(n-24)]
            M[:,24:]*=rng.choice([-1,1],n-24)
            kind='Golay_QR_affine_automorphism__signed_tail_permutation'
        elif i%3==0:
            permutation=rng.permutation(n-24)+24
            M[24:,24:]=np.eye(n-24,dtype=np.int64)[permutation-24]
            signs=1-2*code[int(rng.integers(1,len(code)))]
            M[:24,:24]=np.diag(signs)
            M[:,24:]*=rng.choice([-1,1],n-24)
            kind='Golay_sign_and_signed_tail_permutation'
        else:
            # Two normal-octad coordinates and two tail coordinates: this
            # changes complete midpoint fibres rather than sampled points.
            if extended and i%3==1:
                tail_count=min(4,n-24)
                indices=np.r_[rng.choice(octad,4-tail_count,replace=False),rng.choice(np.arange(24,n),tail_count,replace=False)]
            else:
                indices=np.r_[rng.choice(octad,2,replace=False),rng.choice(np.arange(24,n),2,replace=False)]
            normal=np.zeros(n,np.int64);normal[indices]=rng.choice([-1,1],4)
            M=2*M-np.outer(normal,normal);denominator=2
            kind='tail_triality_reflection' if extended and i%3==1 else 'two_core_two_tail_reflection'
            if i%3==2:
                signs=np.r_[1-2*code[int(rng.integers(1,len(code)))],rng.choice([-1,1],n-24)]
                M*=signs[None,:]
                kind+='__Golay_sign'
        if not np.array_equal(M@M.T,denominator**2*np.eye(n,dtype=np.int64)):
            raise ValueError('Exact isometry')
        yield dict(kind=kind,numerator=M.tolist(),denominator=denominator)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('dimension',type=int)
    ap.add_argument('--trials',type=int,default=12);ap.add_argument('--seed',type=int,default=20261104)
    ap.add_argument('--seconds',type=float,default=600);ap.add_argument('--blockers',type=int,default=8)
    ap.add_argument('--extended',action='store_true')
    ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.dimension not in (27,29,30,31) or args.trials<=0 or args.seconds<=0 or args.blockers<0:
        raise ValueError('Search dimensions/budgets')
    args.out.mkdir(parents=True,exist_ok=False)
    module=load_module(ROOT/'research/discovery_campaign_20261104/competitive28/coupled_search.py')
    paths=list((ROOT/'verifications/improved/dim27-31-joint-refills/data').glob(f'coordinates{args.dimension}_*.npz'))
    if len(paths)!=1:raise ValueError('Unique baseline required')
    source=paths[0];baseline=np.load(source,allow_pickle=False)['points'].astype(np.int16)
    shell=(2*build()).astype(np.int16)
    core={row.tobytes() for row in shell}
    if len(core)!=196560 or not np.all(np.sum(shell.astype(np.int64)**2,axis=1)==128):
        raise ValueError('Regenerated Leech shell')
    start=time.monotonic();records=[];pools=[];deadline=start+args.seconds
    for ordinal,action in enumerate(maps(args.dimension,args.trials,args.seed,args.extended)):
        if time.monotonic()>deadline:break
        M=np.array(action['numerator'],np.int64);d=action['denominator']
        numerator=baseline.astype(np.int64)@M
        if np.any(numerator%d):raise ValueError('Fractional image; this integer search does not round')
        image=(numerator//d).astype(np.int16)
        if not np.all(np.sum(image.astype(np.int64)**2,axis=1)==128):raise ValueError('Image norm')
        final,record,witness=module.exchange(baseline,image,128,blocker_limit=args.blockers,
            max_checks=6000000000,max_edges=5000000,deadline=deadline,core_lookup=core)
        record.update(trial=ordinal,action=action)
        if witness is not None:
            name=f'witness{ordinal:03}.npz';np.savez_compressed(args.out/name,**witness)
            record['witness']=name
            pools.append({'points':witness['restricted_image_points'],'indptr':witness['blocker_indptr'],
                          'blockers':witness['blocker_baseline_indices'],'factor':1})
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in ('action',)}),flush=True)
        (args.out/'report.json').write_text(json.dumps({'dimension':args.dimension,'source':str(source.relative_to(ROOT)),
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline_points':len(baseline),
            'registered_source_geometry':'Inherited for search; fresh complete output replay is mandatory for a gain',
            'records':records,'seed':args.seed,'seconds':time.monotonic()-start,
            'scope':'Finite blocker-restricted whole-isometry images, no global exclusion or novelty claim'},indent=2)+'\n')
        if record.get('gain',0)>0:
            added=image[witness['added_image_indices']]
            proof=module.direct_delta(baseline,final,added)
            np.savez_compressed(args.out/'candidate.npz',points=final,squared_norm=128,
                                inner_product_limit=64,coordinate_denominator=2)
            (args.out/'delta_proof.json').write_text(json.dumps(proof,indent=2)+'\n')
            print('POSITIVE candidate; independent full replay required',flush=True);break
    if pools and not (args.out/'candidate.npz').exists():
        result=module.shared_blocker_search(baseline,pools,max_tests=2000000)
        if result.get('state')=='gain':
            final=result.pop('points');added=result.pop('added')
            np.savez_compressed(args.out/'candidate.npz',points=final,squared_norm=128,
                                inner_product_limit=64,coordinate_denominator=2)
        (args.out/'multi_image.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)

if __name__=='__main__':main()
