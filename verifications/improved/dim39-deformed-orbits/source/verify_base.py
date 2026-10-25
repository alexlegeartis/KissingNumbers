"""Independent exact verifier for the accompanying 39D construction.

Python standard library only. No optimizer, floating point, or saved
numerical distance table is used. Run without -O.
"""
from pathlib import Path
from itertools import combinations,product
import hashlib,json,time
HERE=Path(__file__).resolve().parent
def weight(x):return bin(x).count('1')
def hexrows(name):
    return [int(s,16) for s in (HERE/'source'/name).read_text().splitlines()
            if s and s[0] not in '#$']
def main():
    if not __debug__:raise RuntimeError('Run without -O.')
    start=time.monotonic();cert=json.loads((HERE/'base_certificate.json').read_text())
    assert cert['dimension']==39
    for name,expected in cert['source_sha256'].items():
        assert hashlib.sha256((HERE/'source'/name).read_bytes()).hexdigest()==expected
    generators=hexrows('dense_kernel.txt');reps=hexrows('dense_representatives.txt')
    assert len(generators)==15 and len(reps)==10
    assert all(0<=x<1<<39 for x in generators+reps)
    kernel=[0]
    for g in generators:kernel += [w^g for w in kernel]
    assert len(kernel)==len(set(kernel))==32768
    distance=min(weight(w) for w in kernel if w)
    for a,b in combinations(reps,2):distance=min(distance,min(weight(a^b^w) for w in kernel))
    assert distance>=10
    dense_count=len(reps)*len(kernel)
    print('Dense code:',dense_count,'words; exact minimum distance',distance,flush=True)

    source=(HERE/'source/supports_3424.txt').read_text().splitlines()
    assert len(source)==3424 and all(len(s)==39 and set(s)<={'0','1'} for s in source)
    supports=[tuple(i for i,c in enumerate(s) if c=='1') for s in source]
    for change in cert['support_replacements']:
        i=change['remove_index'];assert list(supports[i])==change['remove']
        supports[i]=tuple(change['add'])
    assert len(set(supports))==3424
    seen=set()
    for s in supports:
        assert len(s)==len(set(s))==8 and all(0<=i<39 for i in s)
        for t in combinations(sorted(s),5):
            assert t not in seen,('Two supports share five coordinates',t)
            seen.add(t)
    print('Updated supports: 3424; every pair intersects in at most four coordinates.',flush=True)

    deleted_pairs=set();added=[]
    for star in cert['heavy_stars']:
        a=star['center'];B=star['leaves']
        assert 0<=a<39 and len(B)==len(set(B))==4 and a not in B
        assert all(0<=b<39 for b in B)
        for b in B:deleted_pairs.add(tuple(sorted((a,b))))
        for center_sign in (-1,1):
            for signs in product((-1,1),repeat=4):
                if sum(s<0 for s in signs)%2:continue
                v=[0]*39;v[a]=4*center_sign
                for b,s in zip(B,signs):v[b]=2*s
                added.append(tuple(v))
    assert len(added)==len(set(added))==16*len(cert['heavy_stars'])
    assert all(sum(x*x for x in v)==32 for v in added)
    assert all(sorted(abs(x) for x in v if x)==[2,2,2,2,4] for v in added)
    max_new=max(sum(x*y for x,y in zip(a,b)) for a,b in combinations(added,2))
    assert max_new<=16
    max_root=-10**9;max_middle=-10**9
    retained_pairs=[p for p in combinations(range(39),2) if p not in deleted_pairs]
    for v in added:
        # Exactly maximize over all four signs of each retained root pair.
        max_root=max(max_root,max(4*(abs(v[i])+abs(v[j])) for i,j in retained_pairs))
        # Exactly maximize over the 128 even sign patterns of each support.
        for s in supports:
            values=[v[i] for i in s]
            score=2*sum(abs(x) for x in values)
            if sum(x<0 for x in values)%2 and all(values):score-=4*min(abs(x) for x in values)
            max_middle=max(max_middle,score)
        # Uniform bound against every +/-1 dense word, using integer squares.
        assert 4*sum(abs(x) for x in v)**2<=32*39
    assert max_root<=16 and max_middle<=16
    assert 4*16**2<=32*39 # Dense versus every middle vector.
    assert 4*8**2<=32*39  # Dense versus every retained root.
    # Middle/middle is certified by the support intersection and even parity.
    # Root/root and middle/root products are at most 16 by their shapes.
    root_count=4*len(retained_pairs);middle_count=128*len(supports)
    total=dense_count+middle_count+root_count+len(added)
    assert total==cert['claimed_size'] and total>cert['previous_size']
    report={'all_checks_passed':True,'dimension':39,'lower_bound':total,
            'dense_words':dense_count,'dense_minimum_distance':distance,
            'middle_vectors':middle_count,'retained_roots':root_count,
            'added_heavy_vectors':len(added),'deleted_root_pairs':sorted(deleted_pairs),
            'maximum_added_added_inner_product':max_new,
            'maximum_added_middle_inner_product':max_middle,
            'maximum_added_root_inner_product':max_root,
            'arithmetic':'Python integers only; no numerical tolerance',
            'certificate_sha256':hashlib.sha256((HERE/'base_certificate.json').read_bytes()).hexdigest(),
            'elapsed_seconds':time.monotonic()-start}
    (HERE/'base_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));print('VERIFIED: tau_39 >= '+str(total)+'.')
if __name__=='__main__':main()
