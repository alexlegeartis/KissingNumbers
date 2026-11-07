#!/usr/bin/env python3
"""Independent exact checker: imports no producer/search/checker modules.
Run with taskset -c 0 env OPENBLAS_NUM_THREADS=1 python independent_audit.py.
Requires Python 3 and NumPy. All products use int16 or wider, never int8.
"""
import argparse, collections, copy, hashlib, itertools, json, os, pathlib, time
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import numpy as np

ROOT=pathlib.Path(__file__).resolve().parent
DEFAULT=ROOT/'certificate32_352832.json'
BASELINE=ROOT/'control32_352800.json'
EXPECTED='69c3f368a4aae5956ca75a8569f761224148831e8a85edd5972a0b3b057af646'

class Rejected(ValueError): pass

def need(condition, reason):
    if not condition: raise Rejected(reason)

def mask(value): return type(value) is int and 0<=value<2**32

def coordinates(s,n):
    return [0 if not ((s>>j)&1) else (-2 if ((n>>j)&1) else 2) for j in range(32)]

def check(d, materialize=True):
    need(d['format']=='signed_octad_kissing_v1','format')
    need(d['dimension']==32,'dimension')
    gen=d['dense_generator']; full=d['full_even_supports']; lines=d['exceptional_antipodal_lines']
    need(len(gen)==17 and all(mask(v) for v in gen),'dense generator format')
    words=[0]
    for g in gen: words += [w^g for w in words]
    need(len(set(words))==2**17,'dense rank')
    wt=collections.Counter(w.bit_count() for w in words)
    need(min(w.bit_count() for w in words if w)>=8,'dense minimum distance')
    need(all(mask(s) and s.bit_count()==8 for s in full),'full support shape')
    need(len(set(full))==len(full),'duplicate full support')
    interfull=collections.Counter()
    for i,s in enumerate(full):
        for t in full[:i]:
            k=(s&t).bit_count(); interfull[k]+=1
            need(k<=4,'full/full intersection')
    ev=[]; masks=[]
    for pair in lines:
        need(isinstance(pair,list) and len(pair)==2,'exceptional pair shape')
        s,n=pair
        need(mask(s) and s.bit_count()==8,'exceptional support shape')
        need(mask(n) and (n&~s)==0,'exceptional negative mask')
        ev += [coordinates(s,n),coordinates(s,s^n)]
        masks += [(s,n),(s,s^n)]
    need(len(set(masks))==len(masks),'duplicate exceptional vector')
    es=sorted({s for s,n in lines}); mix=collections.Counter()
    for s in es:
        for t in full:
            k=(s&t).bit_count(); mix[k]+=1
            need(k<=4,'full/exceptional intersection')
    E=np.asarray(ev,dtype=np.int16).reshape(-1,32)
    need(np.all(np.sum(E*E,axis=1)==32),'exceptional norm')
    gram=E@E.T
    upper=gram[np.triu_indices(len(E),1)]
    need(np.all(upper<=16),'exceptional pair dot')
    hist=dict(zip(*[a.tolist() for a in np.unique(upper,return_counts=True)]))
    patterns=np.asarray([[(-2 if ((b>>j)&1) else 2) for j in range(8)] for b in range(256) if b.bit_count()%2==0],dtype=np.int16)
    pgram=patterns@patterns.T
    need(np.max(pgram[np.triu_indices(128,1)])==16,'full cube internal dot')
    roots=[]
    for a,b in itertools.combinations(range(32),2):
        for x,y in itertools.product((-4,4),repeat=2):
            row=[0]*32; row[a]=x; row[b]=y; roots.append(row)
    R=np.asarray(roots,dtype=np.int16)
    rgram=R@R.T
    need(np.max(rgram[np.triu_indices(len(R),1)])==16,'root pair dot')
    counts={'dense':len(words),'full_even':len(full)*len(patterns),'exceptional':len(E),'roots':len(R)}
    total=sum(counts.values()); need(type(d['points']) is int and d['points']==total,'claimed point count')
    result={'counts':counts,'total':total,'dense_weight_distribution':dict(sorted(wt.items())),
      'dense_rank':17,'dense_minimum_distance':min(k for k in wt if k),
      'full_supports':len(full),'exceptional_supports':len(es),
      'full_full_intersection_histogram':dict(sorted(interfull.items())),
      'full_exceptional_intersection_histogram':dict(sorted(mix.items())),
      'exceptional_unordered_pairs_checked':len(upper),'exceptional_pair_dot_histogram':hist,
      'root_unordered_pairs_checked':len(R)*(len(R)-1)//2,
      'pair_type_maximum_proven_upper_bounds':{'dense_dense':32-2*min(k for k in wt if k),'dense_full':16,'dense_exceptional':16,'dense_root':8,'full_full':16,'full_exceptional':16,'full_root':16,'exceptional_exceptional':int(np.max(upper)),'exceptional_root':16,'root_root':16},
      'global_maximum_distinct_pair_dot':16,'minimum_squared_euclidean_distance':32,'squared_norm':32,
      'pair_coverage':'All ten family pair types covered by exact checks or pointwise shape bounds; not a brute force full Gram matrix.'}
    if materialize:
        D=(1-2*((np.asarray(words,dtype=np.uint32)[:,None]>>np.arange(32,dtype=np.uint32))&1).astype(np.int16))
        F=np.zeros((len(full)*128,32),dtype=np.int16)
        for i,s in enumerate(full):
            idx=[j for j in range(32) if (s>>j)&1]
            F[i*128:(i+1)*128,idx]=patterns
        A=np.concatenate((D,F,E,R))
        need(A.shape==(total,32),'expanded shape')
        norms=np.sum(A*A,axis=1); need(np.all(norms==32),'expanded norms')
        small=A.astype(np.int8)
        unique=len(set(row.tobytes() for row in small))
        need(unique==total,'expanded distinctness')
        result['materialized_coordinate_rows']=len(A)
        result['materialized_distinct_rows']=unique
        result['coordinates_sha256_int8_rowmajor']=hashlib.sha256(small.tobytes()).hexdigest()
    return result

def differences(a,b):
    need(a['dense_generator']==b['dense_generator'],'changed dense generator')
    need(a['full_even_supports']==b['full_even_supports'],'changed full supports')
    def pts(d): return {(s,n) for s,n in d['exceptional_antipodal_lines']} | {(s,s^n) for s,n in d['exceptional_antipodal_lines']}
    old,new=pts(a),pts(b); supports=sorted({s for s,n in old}|{s for s,n in new})
    need({s for s,n in old}=={s for s,n in new},'changed exceptional supports')
    remaining=set(supports); comps=[]
    while remaining:
        todo=[min(remaining)]; comp=set()
        while todo:
            s=todo.pop()
            if s in comp: continue
            comp.add(s); remaining.discard(s)
            todo.extend(t for t in remaining if (s&t).bit_count()>4)
        comps.append({'supports':sorted(comp),'support_count':len(comp),
          'old_points':sum(s in comp for s,n in old),'new_points':sum(s in comp for s,n in new)})
    return {'dense_generator_unchanged':True,'full_supports_unchanged':True,'same_exceptional_supports':len(supports),
      'removed_exceptional_points':len(old-new),'added_exceptional_points':len(new-old),'net_gain':len(new)-len(old),
      'changed_supports':sorted({s for s,n in old^new}),'components':comps}

def negative_controls(original):
    out=[]
    def test(name,edit,expected):
        d=copy.deepcopy(original); edit(d)
        try: check(d,False)
        except Rejected as e:
            need(str(e)==expected,f'control {name}: unexpected rejection {e}')
            out.append({'mutation':name,'rejected':True,'reason':str(e)})
        else: raise Rejected(f'negative control accepted: {name}')
    test('claim one extra point',lambda d:d.__setitem__('points',d['points']+1),'claimed point count')
    test('zero a generator row',lambda d:d['dense_generator'].__setitem__(0,0),'dense rank')
    test('insert weight-one generator',lambda d:d['dense_generator'].__setitem__(0,1),'dense minimum distance')
    test('duplicate full support',lambda d:d['full_even_supports'].__setitem__(0,d['full_even_supports'][1]),'duplicate full support')
    def badfull(d):
        s=d['full_even_supports'][0]; t=(s^(s&-s))|next(1<<i for i in range(32) if not s&(1<<i))
        d['full_even_supports'].append(t); d['points']+=128
    test('append full support intersecting another in seven places',badfull,'full/full intersection')
    test('duplicate exceptional line',lambda d:d['exceptional_antipodal_lines'].__setitem__(0,d['exceptional_antipodal_lines'][1]),'duplicate exceptional vector')
    def outside(d):
        s,n=d['exceptional_antipodal_lines'][0]
        d['exceptional_antipodal_lines'][0][1]=n|next(1<<i for i in range(32) if not s&(1<<i))
    test('put exceptional minus sign outside support',outside,'exceptional negative mask')
    def flip(d):
        s,n=d['exceptional_antipodal_lines'][0];d['exceptional_antipodal_lines'][0][1]=n^(s&-s)
    test('flip one exceptional coordinate',flip,'exceptional pair dot')
    def mixtest(d):
        d['full_even_supports']=[d['exceptional_antipodal_lines'][0][0]]
        d['points']-=128*(len(original['full_even_supports'])-1)
    test('use an exceptional support as a full layer',mixtest,'full/exceptional intersection')
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--certificate',type=pathlib.Path,default=DEFAULT);p.add_argument('--baseline',type=pathlib.Path,default=BASELINE);p.add_argument('--expected-sha256',default=EXPECTED);p.add_argument('--output',type=pathlib.Path,default=pathlib.Path(__file__).with_name('report.json'));args=p.parse_args()
    start=time.monotonic();raw=args.certificate.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    need(sha==args.expected_sha256,'certificate SHA256')
    d=json.loads(raw);basebytes=args.baseline.read_bytes();base=json.loads(basebytes)
    report={'certificate':str(args.certificate.relative_to(ROOT)),'certificate_sha256':sha,
      'baseline':str(args.baseline.relative_to(ROOT)),'baseline_sha256':hashlib.sha256(basebytes).hexdigest(),
      'certificate_audit':check(d),'baseline_audit':check(base),'delta':differences(base,d),
      'negative_controls':negative_controls(d),'cpu_affinity':sorted(os.sched_getaffinity(0)),
      'python_numpy_version':np.__version__,'independence':'New standalone implementation; no producer checker/search imports or calls.'}
    report['status']='PASS';report['elapsed_seconds']=time.monotonic()-start
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'total':report['certificate_audit']['total'],'delta':report['delta']['net_gain'],'negative_controls':len(report['negative_controls']),'seconds':report['elapsed_seconds'],'report':str(args.output)},indent=2))
if __name__=='__main__': main()
