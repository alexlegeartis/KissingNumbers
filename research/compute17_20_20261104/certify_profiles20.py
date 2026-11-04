#!/usr/bin/env python3
"""Exact 19448 ceiling for the saved 267-profile, 4096-word finite family.

No optimizer and no floating point enter this proof. This is NOT a bound on
arbitrary nonuniform profiles or on the unrestricted kissing number.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import numpy as np
from orbits17 import DATA
from prepare20 import basis_of
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
from golay import golay24


def main():
    source=DATA/'profiles20.json'
    profiles=[tuple(map(int,p)) for p in json.loads(source.read_text())['profiles']]
    assert len(profiles)==len(set(profiles))
    assert all(len(p)==20 and min(p)>=0 and math.gcd(*p)==1 for p in profiles)
    G=golay24().astype(int)
    words=[sum(int(x)<<i for i,x in enumerate(row)) for row in G[:,4:]]
    assert len(set(words))==4096
    blocks=[tuple(np.flatnonzero(row[4:])) for row in G if sum(row[:4])==0 and sum(row)==8]
    assert len(blocks)==130
    norms=[sum(x*x for x in p) for p in profiles]
    # All 267*4096 vectors are compatible with the 17400 fixed root/octad points.
    # v = sqrt(8/Q)*p*(-1)^c. Since c is orthogonal to each shortened octad,
    # its signs disagree with every odd octad pattern at least once.
    for p,Q in zip(profiles,norms):
        assert all(2*(p[i]+p[j])**2<=Q for i,j in itertools.combinations(range(20),2))
        assert all((sum(p[i] for i in block)-2*min(p[i] for i in block))**2<=2*Q for block in blocks)
    critical=sorted(w for w in words if w.bit_count()==4)
    assert len(critical)==5
    assert all((a&b)==0 for a,b in itertools.combinations(critical,2))
    assert sum(critical)==(1<<20)-1
    # At most one point can use a given sign word, regardless of its profile.
    for i,p in enumerate(profiles):
        for j,q in enumerate(profiles):
            dot=sum(a*b for a,b in zip(p,q))
            assert dot>0 and 4*dot*dot>norms[i]*norms[j]
    # If a word's profile permits a selected neighbor at critical difference s,
    # mark that direction. Every saved profile permits at most one direction,
    # even when its neighbor can use any of the other 266 profiles.
    enabled=[]
    for i,p in enumerate(profiles):
        directions=[]
        for color,word in enumerate(critical):
            for j,q in enumerate(profiles):
                dot=sum(a*b*(-1 if word>>k&1 else 1) for k,(a,b) in enumerate(zip(p,q)))
                if dot<=0 or 4*dot*dot<=norms[i]*norms[j]:
                    directions.append(color)
                    break
        assert len(directions)<=1
        enabled.append(directions)
    # Partition the 4096 words into 1024 squares using two critical directions.
    # Three selected words in a square would give its middle word two selected
    # critical neighbors, impossible by the preceding exact profile checks.
    a,b=critical[:2];unseen=set(words);squares=[]
    while unseen:
        w=min(unseen);square={w,w^a,w^b,w^a^b}
        assert len(square)==4 and square<=unseen
        unseen-=square;squares.append(sorted(square))
    assert len(squares)==1024
    upper=2*len(squares)
    # Matching lower bound: the original uniform profile and a linear
    # 2048-word subcode with minimum Hamming distance five.
    assert tuple([1]*20) in profiles
    F={w for w,row in zip(words,G) if row[0]==0}
    assert len(F)==2048 and 0 in F
    fb=basis_of(F);span={0}
    for gen in fb:span|={w^gen for w in span}
    assert span==F and len(fb)==11
    assert min(w.bit_count() for w in F if w)==5
    assert len(F)==upper
    report=dict(status='PASS',method='integer profile comparisons and a partition into squares',
                profiles=len(profiles),sign_words=4096,dense_candidate_points=len(profiles)*4096,
                fixed_base=17400,max_added_layer=upper,exact_family_ceiling=17400+upper,
                square_count=len(squares),critical_word_masks=critical,
                profiles_without_enabled_critical_direction=sum(not x for x in enabled),
                profiles_with_one_enabled_critical_direction=sum(bool(x) for x in enabled),
                enabled_critical_directions=enabled,
                profile_file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                scope='Only the saved 267 profiles on the 4096 punctured Golay sign words with the complete fixed 17400-point base.')
    (DATA/'profile_ceiling_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='enabled_critical_directions'}),flush=True)


if __name__=='__main__':main()
