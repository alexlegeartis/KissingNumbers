#!/usr/bin/env python3
"""Independent audit of the pointwise-star 370408 continuation.

Imports only this auditor's independent_core.py, never a producer/search/proof module.
Root deletions and stars are INDIVIDUAL signed rows, with no antipodal closure.
"""
import argparse
import copy
import hashlib
import gzip
import itertools
import json
from pathlib import Path
import time

import independent_core as own
import numpy as np

PATH = ''
PINNED = {
    'certificate33_370408_explicit.json': '5686dfcd2b891bf6d9be6b01013d8205774aa1a19c02a2c45bdfd86194314a60',
    'delta_370408.json': 'bd494930e97ebfa8d72c73988b4f67a811deeb8ea998298a651f0df698c17f9c',
}
req = own.require


def point_schema(d):
    own.schema(d)
    req(d['format'] == 'signed_octad_star_points_explicit_v1', 'wrong pointwise format')
    req(d['dimension'] == 33, 'wrong ambient dimension')
    req('stars' not in d and 'removed_root_lines' not in d, 'mixed line/point star/root semantics')
    for field, shape in [('removed_root_signed_points', [4,4]), ('star_signed_points', [2,2,2,2,4])]:
        values = d[field]
        req(len(values) == len(set(map(tuple, values))), 'duplicate individual row in '+field)
        for x in values:
            req(len(x) == 33 and all(type(y) is int for y in x), 'bad integer row length in '+field)
            req(sorted(abs(y) for y in x if y) == shape, 'bad magnitude shape in '+field)
            req(sum(y*y for y in x) == 32, 'bad norm in '+field)


def desc(x, kind):
    s = sum(1 << i for i, v in enumerate(x) if v)
    neg = sum(1 << i for i, v in enumerate(x) if v < 0)
    center = next(i for i, v in enumerate(x) if abs(v) == 4) if kind == 2 else 0
    return kind, s, neg, center


def reconstruct(d):
    point_schema(d)
    core = copy.deepcopy(d)
    del core['removed_root_signed_points']
    del core['star_signed_points']
    core['points'] = 131072 + 128*len(core['full_even_supports']) + len(core['exceptional_signed_points']) + 2112
    rec = own.reconstruct(core)
    removal = {desc(x, 1) for x in d['removed_root_signed_points']}
    universe = set(rec['descriptors']['root'])
    req(removal <= universe, 'removed root not in complete root universe')
    for x in removal:
        del rec['sparse'][x]
    rec['descriptors']['root'] = [x for x in rec['descriptors']['root'] if x not in removal]
    rec['arrays']['root'] = own.rows(rec['descriptors']['root'])
    rec['descriptors']['star'] = []
    for i, x in enumerate(d['star_signed_points']):
        dd = desc(x, 2)
        req(dd not in rec['sparse'], 'duplicate reconstructed pointwise star')
        req(np.array_equal(own.row(dd), x), 'star coordinate reconstruction differs')
        rec['sparse'][dd] = ['star_point', i]
        rec['descriptors']['star'].append(dd)
    rec['arrays']['star'] = np.asarray(d['star_signed_points'], dtype=np.int8).reshape((-1,33))
    rec['rows'] = np.concatenate(list(rec['arrays'].values()))
    rec['rowset'] = {x.tobytes() for x in rec['rows']}
    req(len(rec['rowset']) == len(rec['rows']), 'global duplicate pointwise row')
    req(len(rec['rows']) == d['points'], 'claimed pointwise count differs from rows')
    req(np.all(np.sum(rec['rows'].astype(np.int64)**2, axis=1) == 32), 'pointwise row norm differs')
    rec['data'] = d
    return rec


def setrows(a):
    return {tuple(int(v) for v in x) for x in a}


def delta_proof(c, base, actual, delta, pool):
    ca, ba, aa = c['arrays'], base['arrays'], actual['arrays']
    for cat in ('dense','full','explicit'):
        req(np.array_equal(ca[cat], ba[cat]), 'strong candidate changes '+cat)
    cs, bs = setrows(ca['star']), setrows(ba['star'])
    cr, br = setrows(ca['root']), setrows(ba['root'])
    changes = {'added_star_points': cs-bs, 'removed_star_points': bs-cs,
               'restored_root_points': cr-br, 'newly_removed_root_points': br-cr}
    for key, value in changes.items():
        req(len(delta[key]) == len(set(map(tuple, delta[key]))), 'duplicate delta row in '+key)
        req(set(map(tuple, delta[key])) == value, 'exact delta differs for '+key)
    req(delta['source_certificate_sha256'] == own.FILES['candidate'][1], 'delta source hash differs')
    req(delta['new_certificate_sha256'] == PINNED['certificate33_370408_explicit.json'], 'delta candidate hash differs')
    req({k:len(v) for k,v in changes.items()} == {'added_star_points':148,'removed_star_points':32,'restored_root_points':0,'newly_removed_root_points':100}, 'wrong pointwise delta counts')
    req(c['rowset'] == (base['rowset']-{np.asarray(x,dtype=np.int8).tobytes() for x in (bs-cs)|(br-cr)}) | {np.asarray(x,dtype=np.int8).tobytes() for x in (cs-bs)|(cr-br)}, 'full coordinate set delta does not reproduce candidate')
    req(pool['source_sha256'] == own.FILES['candidate'][1], 'pool source hash differs')
    independent_roots = []
    for i,j in itertools.combinations(range(33),2):
        for a,b in itertools.product((-1,1),repeat=2):
            x = [0]*33
            x[i],x[j] = 4*a,4*b
            independent_roots.append(x)
    req(pool['all_roots'] == independent_roots, 'pool root ID ordering differs from reconstructed universe')
    for key, values, universe in [('new_star_ids',c['data']['star_signed_points'],pool['star_vectors']),
                                  ('deleted_root_ids',c['data']['removed_root_signed_points'],independent_roots)]:
        ids = delta[key]
        req(len(ids) == len(set(ids)), 'duplicate pool selection ID')
        req(all(type(i) is int and 0 <= i < len(universe) for i in ids), 'bad pool selection ID')
        req([universe[i] for i in ids] == values, key+' does not reproduce certificate')
    unpaired = {cat:sum(tuple(-int(y) for y in x) not in setrows(rows) for x in rows)
                for cat, rows in ca.items() if cat in ('star','root','explicit')}
    req(unpaired == {'explicit':372,'root':12,'star':116}, 'pointwise missing-antipode counts differ')
    req(len(c['rowset']-actual['rowset']) == 1796 and len(actual['rowset']-c['rowset']) == 1764, 'actual old control coordinate delta differs')
    return {'changes_vs_370392': {k:len(v) for k,v in changes.items()},
            'gain_vs_370392':len(c['rows'])-len(base['rows']),
            'gain_vs_actual_370376':len(c['rows'])-len(actual['rows']),
            'actual_old_control_added':len(c['rowset']-actual['rowset']),
            'actual_old_control_removed':len(actual['rowset']-c['rowset']),
            'unchanged_dense':len(ca['dense']), 'unchanged_full':len(ca['full']),
            'unchanged_exceptional':len(ca['explicit']), 'all_roots_before_deletion':2112,
            'individual_deleted_roots':156, 'retained_roots':len(ca['root']),
            'individual_stars':len(ca['star']), 'missing_antipodes_by_category':unpaired,
            'pool_source_id_checks':len(delta['new_star_ids'])+len(delta['deleted_root_ids']),
            'scope':'Pool IDs are provenance checks only. No completeness, solver, MILP optimum, or search-graph result is trusted.'}


def controls(d, rec, raw):
    outcomes = []
    def reject(name, fn, expected=None):
        try:
            fn()
        except ValueError as e:
            req(expected is None or expected in str(e), name+' rejected for unintended reason: '+str(e))
            outcomes.append({'name':name,'rejected':True,'reason':str(e)})
            return
        raise ValueError(name+': corruption was accepted')
    def mutation(fn):
        x=copy.deepcopy(d)
        fn(x)
        return x
    reject('corrupt stronger candidate byte hash',lambda:own.checked_json(raw+b' ', PINNED['certificate33_370408_explicit.json']),'SHA-256')
    x=mutation(lambda x:x['star_signed_points'].append(x['star_signed_points'][0]))
    reject('duplicate individual star',lambda:point_schema(x),'duplicate individual')
    x=mutation(lambda x:x['removed_root_signed_points'].append(x['removed_root_signed_points'][0]))
    reject('duplicate individual root deletion',lambda:point_schema(x),'duplicate individual')
    x=mutation(lambda x:x['star_signed_points'][0].pop())
    reject('bad individual star dimension',lambda:point_schema(x),'row length')
    x=mutation(lambda x:x['star_signed_points'][0].__setitem__(0,-3))
    reject('bad individual star magnitude',lambda:point_schema(x),'magnitude shape')
    x=mutation(lambda x:x['removed_root_signed_points'][0].__setitem__(0,-2))
    reject('bad individual root magnitude',lambda:point_schema(x),'magnitude shape')
    x=mutation(lambda x:x.__setitem__('removed_root_lines',[]))
    reject('mixed root line and signed-point schema',lambda:point_schema(x),'mixed line/point')
    x=mutation(lambda x:x.__setitem__('points',x['points']+1))
    reject('incorrect stronger total',lambda:reconstruct(x),'claimed pointwise count')
    x=mutation(lambda x:x['dense_generator'].__setitem__(0,x['dense_generator'][1]))
    reject('deficient dense code rank',lambda:own.dense(x),'deficient rank')
    x=mutation(lambda x:x['dense_coordinates'].append(32))
    reject('wrong dense support size33',lambda:point_schema(x),'exactly 32')
    x=mutation(lambda x:x['exceptional_signed_points'][0].__setitem__(1,1<<32))
    reject('exceptional mask outside support',lambda:point_schema(x),'outside octad')
    # Adding antipodes of deletion rows deletes 12 additional retained roots.
    # It changes the construction and must not silently preserve the count.
    removal=set(map(tuple,d['removed_root_signed_points']))
    completed=removal|{tuple(-y for y in x) for x in removal}
    req(len(completed)==168,'expected individual root deletion asymmetry absent')
    x=mutation(lambda x:x.__setitem__('removed_root_signed_points',[list(v) for v in sorted(completed)]))
    reject('incorrect antipodal completion of root deletions',lambda:reconstruct(x),'claimed pointwise count')
    restored=np.asarray(d['removed_root_signed_points'][0],dtype=np.int8)
    reject('restore one forbidden individual root',lambda:own.pair_products(restored[None,:],rec['arrays']['star'],'restored root/star'),'24')
    witness={}
    for cat in ('star','root'):
        vals=rec['arrays'][cat]
        present={x.tobytes() for x in vals}
        chosen=None
        for i,x in enumerate(vals):
            if (-x).tobytes() in present:
                continue
            for other in ('explicit','root','star'):
                dots=(-x).astype(np.int64)@rec['arrays'][other].astype(np.int64).T
                if len(dots) and int(dots.max())>16:
                    j=int(dots.argmax())
                    chosen={'source_category':cat,'source_index':i,'other_category':other,'other_index':j,
                            'source_row':[int(v) for v in x], 'other_row':[int(v) for v in rec['arrays'][other][j]],
                            'after_antipode_product':int(dots[j]), 'before_product':-int(dots[j])}
                    reject('complete missing '+cat+' antipode',lambda x=x,other=other:own.pair_products((-x)[None,:],rec['arrays'][other],'completed '+cat+'/'+other),'24')
                    break
            if chosen:
                break
        req(chosen is not None,'missing antipode conflict witness absent for '+cat)
        witness[cat]=chosen
    # The octad asymmetry remains part of the stronger construction too.
    ep=rec['arrays']['explicit']
    present={x.tobytes() for x in ep}
    chosen=None
    for i,x in enumerate(ep):
        if (-x).tobytes() in present:
            continue
        dots=(-x).astype(np.int64)@ep.astype(np.int64).T
        if int(dots.max())>16:
            j=int(dots.argmax())
            chosen={'source_index':i,'other_index':j,'after_antipode_product':int(dots[j]),'before_product':-int(dots[j])}
            reject('complete missing exceptional antipode',lambda x=x:own.pair_products((-x)[None,:],ep,'completed exceptional'),'24')
            break
    req(chosen is not None,'missing exceptional antipode conflict witness absent')
    witness['explicit']=chosen
    return {'rejection_count':len(outcomes),'results':outcomes,'missing_antipode_conflicts':witness,
            'semantic_note':'All individual signed rows are interpreted literally. Existing -24 products are allowed; their missing antipodes produce forbidden +24 products.'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--workspace',type=Path,default=own.HERE)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    loaded={k:own.checked_json((args.workspace/PATH/k).read_bytes(),v) for k,v in PINNED.items()}
    d=loaded['certificate33_370408_explicit.json']
    req(d['points']==370408,'wrong target count')
    req(d['dense_coordinates']==list(range(32)) and d['dense_coset_shift']==0,'wrong fixed dense placement')
    req(d['construction']['baseline_sha256']==own.FILES['candidate'][1],'wrong construction baseline hash')
    rec=reconstruct(d)
    proof=own.geometry(rec,'stronger370408')
    req(proof['counts']=={'dense':131072,'full':235136,'explicit':2064,'root':1956,'star':180},'wrong stronger category counts')
    req(proof['rank_mod101']==33 and proof['missing_antipodes']==500,'rank or antipodal count differs')
    sources={k:own.checked_json((args.workspace/own.LOCAL_FILES[k]).read_bytes(),own.FILES[k][1]) for k in ('candidate','classic','old_trade')}
    base=own.reconstruct(sources['candidate'])
    actual=copy.deepcopy(sources['classic'])
    actual.update({k:sources['old_trade'][k] for k in ('points','stars','removed_root_lines')})
    actual=own.reconstruct(actual)
    poolbytes=gzip.decompress((args.workspace/PATH/'pool.json.gz').read_bytes())
    req(hashlib.sha256(poolbytes).hexdigest()=='7827321d836a37965b1d3d0dd884cf7d5b8fb26532d37b1c63a4dc4dbbbf3fa9', 'decompressed provenance pool SHA-256 mismatch')
    delta=delta_proof(rec,base,actual,loaded['delta_370408.json'],json.loads(poolbytes))
    result={'status':'PASS','dimension':33,'lower_bound':370408,'candidate':proof,'delta':delta,
            'hashes':PINNED,'pool_sha256':hashlib.sha256(poolbytes).hexdigest(),
            'independence':'Only this independent auditor module is imported. No producer proof, solver, graph, matching, or search function is imported or called.',
            'scope':'Direct finite geometry and exact deltas: +16 versus370392; +32 versus the actual370376 repository control. No worldwide-priority or optimality claim.'}
    (args.output/'geometry_pass.json').write_text(json.dumps(result,indent=2)+'\n')
    print('MATHEMATICAL PASS: 370408 distinct norm-32 rows in R^33; ordinary products <=16; +32 over actual old control.',flush=True)
    result['controls']=controls(d,rec,(args.workspace/PATH/'certificate33_370408_explicit.json').read_bytes())
    result['elapsed_seconds']=time.monotonic()-start
    result['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['independent_core_sha256']=hashlib.sha256((own.HERE/'independent_core.py').read_bytes()).hexdigest()
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('AUDIT PASS including '+str(result['controls']['rejection_count'])+' controls; '+str(args.output/'report.json'),flush=True)


if __name__=='__main__':
    main()
