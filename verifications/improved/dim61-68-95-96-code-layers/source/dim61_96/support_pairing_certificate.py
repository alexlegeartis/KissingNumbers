#!/usr/bin/env python3
"""Independently verify optimized QR48 doubling and its best two-row section."""
import hashlib
import itertools
import json
import time
import numpy as np
from support_doubling import qr48, expand
from section95_search import best_rs_section
from weighted_layers import HERE, sign24, kernel, matrix, packed, all_binary_words, popcounts64


def main():
    start=time.perf_counter()
    values=list(map(int,(HERE/'support_pairing_best.txt').read_text().split()))
    claimed,permutation=values[0],values[1:]
    assert sorted(permutation)==list(range(48))
    generator,words,histogram=qr48()  # Re-enumerates every one of the 2^24 words.
    remapped=np.zeros(len(words),dtype=np.uint64)
    for new,old in enumerate(permutation):
        remapped |= ((words>>np.uint64(old))&np.uint64(1))<<np.uint64(new)
    good=remapped[(remapped & (remapped>>np.uint64(1)) & np.uint64(0x555555555555))==0]
    assert len(good)==claimed
    additions=[expand(int(w),0) for w in remapped]
    additions += [expand(int(w),p) for p in (1,2) for w in good]
    assert len(additions)==len(set(additions))==17296+2*len(good)
    assert all(bin(w).count('1')==24 for w in additions)
    incidence=np.array([[(w>>c)&1 for c in range(96)] for w in additions],dtype=np.int16)
    row_counts=incidence.reshape(-1,24,4).sum(axis=2)
    assert np.all(np.count_nonzero(row_counts,axis=1)<=12)
    assert np.all((row_counts[17296:]==0)|(row_counts[17296:]==2))
    # Same-partition intersections follow from the exhaustively checked QR48 distance.
    # Between partitions each occupied cell pair meets an unlike pair in one cell;
    # one side has no double row, so the sum is at most 12.
    sign_G,_=sign24()
    dual=all_binary_words(matrix(kernel(packed(sign_G),24),24))
    weights=popcounts64(dual);assert min(weights[1:])==8
    dual8=int(dual[np.flatnonzero(weights==8)[0]])
    signs=all_binary_words(sign_G)
    retained=np.array([0 if j%2 else np.count_nonzero(popcounts64(signs & np.uint64(dual8 if j==8 else (1<<j)-1))==j//2)
                       for j in range(9)],dtype=np.int64)
    best=(-1,None,None)
    for ra,rb in itertools.combinations(range(24),2):
        hist=np.bincount(row_counts[:,ra]+row_counts[:,rb],minlength=9)
        count=int(hist@retained)
        if count>best[0]:best=(count,(ra,rb),hist)
    rs=best_rs_section()
    terms95={'dense_section':4697620480,'transversal24_section':2147483648,
             'new_support24_section':best[0],'support6_section':rs['points'],'ordinary_poles':176}
    terms96={'dense':1<<33,'transversal24':1<<32,'new_support24':len(additions)*16384,
             'support6':1<<21,'poles':192}
    coordinates=[4*r+c for r in best[1] for c in range(4)]
    # Check the selected section against the actual dense generator, not a presumed symmetry.
    from weighted_layers import dense96,quaternary24,rref
    q4,_=quaternary24();dense=dense96(sign_G,q4)
    projection=rref([sum(int(row[c])<<j for j,c in enumerate(coordinates)) for row in dense],8)[0]
    projected=all_binary_words(matrix(projection,8))
    assert len(projection)==7 and np.all(popcounts64(projected)%2==0)
    assert np.count_nonzero(popcounts64(projected)==4)*(1<<(33-7))==terms95['dense_section']
    path=HERE/'support_pairing_certificate.npz'
    np.savez_compressed(path,qr48_generator=generator,permutation=np.array(permutation,dtype=np.uint8),
        remapped_words=remapped,no_double_words=good,
        added_supports_low=np.array([w&((1<<64)-1) for w in additions],dtype=np.uint64),
        added_supports_high=np.array([w>>64 for w in additions],dtype=np.uint64),
        section_coordinates=np.array(coordinates,dtype=np.uint8),sign_generator=sign_G)
    report={'status':'Verified: exhaustive QR48 distance and exact support/section counts',
        'no_double_supports':len(good),'added_supports':len(additions),
        'dimensions':{'95':{'bound':sum(terms95.values()),'previous_bound':7077231372,
            'gain':sum(terms95.values())-7077231372,'terms':terms95,'section_coordinates':coordinates,
            'support_intersection_histogram':best[2].tolist(),'retained_sign_counts':retained.tolist()},
            '96':{'bound':sum(terms96.values()),'previous_bound':13264027840,
            'gain':sum(terms96.values())-13264027840,'terms':terms96}},
        'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'QR48_histogram':{str(i):int(c) for i,c in enumerate(histogram) if c},
        'seconds':time.perf_counter()-start,
        'scope':'Constructive improvements over the preceding repository constructions; no literature-priority claim.',
        'fast_check':'Verify the 48-coordinate code once, then use occupancy/intersection inequalities and exact sign projection counts; no enumeration of billions of kissing points.',
        'uniform_partition_obstruction':'For every QR48 word containing a double row, its copies in distinct uniform partitions have support intersection >12. Thus the three uniform partitions alone allow at most 17296+2*N_no_double supports.'}
    (HERE/'support_pairing_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
