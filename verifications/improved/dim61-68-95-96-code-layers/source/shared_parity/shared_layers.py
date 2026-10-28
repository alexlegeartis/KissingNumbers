"""Shared parity repairs, followed by exact transfers to dimensions 61--68.

Every sparse point has 16 nonzero entries +/-2 and squared norm 64.
The archives specify extended Hamming duals and affine shifts, not search
scores.  Intersections are checked in packed and coordinate representations.
"""
from fractions import Fraction as Q
import hashlib
import itertools
import json
from math import comb
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'dim61_96'))
from weighted_layers import (build_code64_rows, collect_min_words, impose_zeros,
                             kernel, matrix, packed, popcounts64, rref, sqs64)
from binary_62_code import DUAL_SUPPORT, EXPONENTS, keep_coordinates
sys.path.insert(0, str(HERE.parent/'significant49_72'))
from compressed_layers import hex_geometry
from latitude_compression import geometry, quad_leq

bits = getattr(int, 'bit_count', lambda word: bin(word).count('1'))
FOUR = (24, 42, 53, 55)
PREVIOUS = {61:87056170, 62:120620044, 63:182537988, 64:332170112,
            65:332238338, 66:338302340, 67:351420806, 68:371552152}


def words(rows):
    out = [0]
    for row in rows:
        out += [word ^ row for word in out]
    return out


def coordinates(support):
    return [c for c in range(64) if support >> c & 1]


def local(support, mask):
    return sum((mask >> c & 1) << i for i, c in enumerate(coordinates(support)))


def expand(support, mask):
    return sum((mask >> i & 1) << c for i, c in enumerate(coordinates(support)))


def solve(equations):
    """Exact sparse binary elimination, returning one solution or None."""
    pivots = {}
    for row, rhs in equations:
        while row:
            pivot = (row & -row).bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = (row, rhs)
                break
            other, value = pivots[pivot]
            row, rhs = row ^ other, rhs ^ value
        if not row and rhs:
            return None, len(pivots)
    value = 0
    for pivot, (row, rhs) in sorted(pivots.items(), reverse=True):
        if bits(row & value) % 2 != rhs:
            value ^= 1 << pivot
    assert all(bits(row & value) % 2 == rhs for row, rhs in equations)
    return value, len(pivots)


def dual_data(support, checks, shift):
    assert checks[0] == 65535 and len(checks) == len(rref(checks,16)[0]) == 5
    assert 0 <= shift < 2**16
    dual = words(checks)
    assert len(set(dual)) == 32 and {bits(h) for h in dual} == {0,8,16}
    # The 16 distinct odd syndrome columns prove no <=3 columns sum to zero.
    columns = [sum((check >> c & 1) << j for j, check in enumerate(checks))
               for c in range(16)]
    assert len(set(columns)) == 16 and all(column & 1 for column in columns)
    code = words(kernel(checks,16))
    assert len(code) == 2048 and min(bits(word) for word in code if word) == 4
    return dict(checks=checks, shift=shift, code=code,
                span={expand(support,h): i for i,h in enumerate(dual)},
                phase=expand(support,shift))


def load_family(name, old, generator):
    path = HERE/(name+'.npz')
    with np.load(path, allow_pickle=False) as archive:
        z = {key: archive[key].copy() for key in archive.files}
    assert set(z) == {kind+'_'+field for kind in ('new','old')
                      for field in ('supports','duals','shifts')}
    assert all(value.dtype.kind == 'u' for value in z.values())
    new = list(map(int,z['new_supports']))
    modified = list(map(int,z['old_supports']))
    assert len(new) == len(set(new)) and len(modified) == len(set(modified))
    assert set(modified) <= set(old) and not set(new).intersection(old)
    assert all(bits(s) == 16 for s in new+modified)
    assert all(bits(row & s) % 2 == 0 for row in generator for s in new)
    data = {}
    for kind, supports in (('new',new),('old',modified)):
        assert z[kind+'_duals'].shape == (len(supports),5)
        assert z[kind+'_shifts'].shape == (len(supports),)
        for s, d, phase in zip(supports,z[kind+'_duals'],z[kind+'_shifts']):
            data[s] = dual_data(s,list(map(int,d)),int(phase))
    old_array = np.array(old,dtype=np.uint64)
    old_coords = np.array([[(s>>c)&1 for c in range(64)] for s in old],dtype=np.float32)
    new_array = np.array(new,dtype=np.uint64)
    new_coords = np.array([[(s>>c)&1 for c in range(64)] for s in new],dtype=np.float32)
    conflicts, touched, new_pairs = [], set(), 0
    for i,a in enumerate(new):
        packed_counts = popcounts64(old_array & np.uint64(a))
        # All operations are exact: sums of <=64 binary entries, below 2^24.
        physical_counts = old_coords @ new_coords[i]
        assert np.array_equal(packed_counts,physical_counts)
        bad = np.flatnonzero(packed_counts > 8)
        assert len(bad) == 8 and int(packed_counts.max()) == 10
        for j in bad:
            b = old[int(j)]
            assert b in data
            overlap = a & b
            shared = [h for h in data[a]['span'] if bits(h) == 8
                      and h & overlap == h and h in data[b]['span']]
            assert len(shared) == 1
            h = shared[0]
            assert bits(data[a]['phase'] & h) % 2 != bits(data[b]['phase'] & h) % 2
            assert 4*(10-2) == 32
            conflicts.append((a,b,h))
            touched.add(b)
        if i+1 < len(new):
            counts = popcounts64(new_array[i+1:] & np.uint64(a))
            assert np.array_equal(counts,new_coords[i+1:] @ new_coords[i])
            assert int(counts.max()) <= 8
            new_pairs += len(counts)
    assert touched == set(modified)
    assert new_pairs == len(new)*(len(new)-1)//2
    receipt = dict(new_supports=len(new), modified_old_supports=len(modified),
        old_new_pairs_checked=len(old)*len(new), new_pairs_checked=new_pairs,
        parity_repairs=len(conflicts), certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return dict(new=new,data=data,conflicts=conflicts,receipt=receipt)


def odd_phases(family, removed):
    """Re-solve ALL surviving parity repairs with odd total signs, for 61D."""
    supports = sorted(s for s in family['data'] if not(s & removed))
    position = {s:i for i,s in enumerate(supports)}
    equations, constraints = [], []
    for a,b,h in family['conflicts']:
        if a & removed or b & removed:
            continue
        ca,cb = family['data'][a]['span'][h],family['data'][b]['span'][h]
        row = ((ca>>1) << (4*position[a])) ^ ((cb>>1) << (4*position[b]))
        # Syndrome bit zero is fixed to one, since the first check is all ones.
        rhs = 1 ^ (ca&1) ^ (cb&1)
        equations.append((row,rhs))
        constraints.append((a,b,h))
    value,rank = solve(equations)
    assert value is not None
    shifts = {}
    for s,i in position.items():
        target = 1 | (((value >> (4*i)) & 15) << 1)
        checks = family['data'][s]['checks']
        equations16 = [(check,target>>j&1) for j,check in enumerate(checks)]
        shift,_ = solve(equations16)
        assert shift is not None and bits(shift)%2 == 1
        assert all(bits(shift & check)%2 == (target>>j&1) for j,check in enumerate(checks))
        shifts[s] = expand(s,shift)
    assert all(bits(shifts[a]&h)%2 != bits(shifts[b]&h)%2 for a,b,h in constraints)
    return dict(supports=len(supports), equations=len(equations), rank=rank)


def best_four(supports):
    masks = []
    for c in range(64):
        bitmap = bytearray((len(supports)+7)//8)
        for j,s in enumerate(supports):
            if not(s>>c&1): bitmap[j//8] |= 1 << (j%8)
        masks.append(int.from_bytes(bitmap,'little'))
    best,choice = -1,None
    for a,b,c,d in itertools.combinations(range(64),4):
        count = bits(masks[a]&masks[b]&masks[c]&masks[d])
        if count > best: best,choice = count,(a,b,c,d)
    assert choice == FOUR and best == 10261
    return dict(choices=comb(64,4), maximum=best, selected_coordinates=list(choice))


def section_count(s,h,data):
    mask,m = s & h,bits(s & h)
    if m % 2: return 0
    if m < 8: return comb(m,m//2)*(1 << (11-m))
    assert m == 8
    if s not in data: return 1120  # Freely put this octad first on an unchanged support.
    if mask in data[s]['span']:
        return 0 if bits(data[s]['phase'] & mask)%2 else 1120
    return 560


def check_section_counts(family, old, section):
    # Independent enumeration, rather than the projection-rank formula above.
    standard = [65535]+[sum((i>>j&1)<<i for i in range(16)) for j in range(4)]
    standard_words = words(kernel(standard,16))
    for s in old+family['new']:
        m = bits(s & section)
        if s in family['data']:
            d = family['data'][s]
            mask = local(s,section)
            count = sum(bits((word ^ d['shift']) & mask) == m//2 for word in d['code']) if m%2 == 0 else 0
        else:
            mask = (1<<m)-1
            count = sum(bits(word & mask) == m//2 for word in standard_words) if m%2 == 0 else 0
        assert count == section_count(s,section,family['data'])


def negative_control(family):
    a,b,h = family['conflicts'][0]
    A,B = family['data'][a],family['data'][b]
    ca,cb = coordinates(a),coordinates(b)
    equations = [(check,bits(check&A['shift'])%2) for check in A['checks']]
    equations += [(check<<16,bits(check&B['shift'])%2) for check in B['checks']]
    equations += [((1<<ca.index(c)) | (1<<(16+cb.index(c))),0) for c in coordinates(a&b)]
    assert solve(equations)[0] is None
    corrupted = B['shift'] ^ (1 << cb.index((h & -h).bit_length()-1))
    bad = [(check,bits(check&A['shift'])%2) for check in A['checks']]
    bad += [(check<<16,bits(check&corrupted)%2) for check in B['checks']]
    bad += equations[10:]
    witness,_ = solve(bad)
    assert witness is not None
    left,right = witness & 65535,witness>>16
    ip = sum(4*(1-2*(left>>ca.index(c)&1))*(1-2*(right>>cb.index(c)&1))
             for c in coordinates(a&b))
    assert ip == 40 > 32
    return dict(rejected_original_matching_system=True, corrupted_phase_inner_product=ip)


def main():
    if not __debug__: raise SystemExit('Run without -O: assertions are required.')
    generator = build_code64_rows()
    raw,histogram = collect_min_words(matrix(generator,64),16)
    old = sorted(map(int,raw))
    assert len(old) == len(set(old)) == 30828
    assert all(bits(a&b)%2 == 0 for a in generator for b in generator)
    with np.load(HERE.parent/'dim61_96/code64_supports.npz',allow_pickle=False) as archive:
        assert old == sorted(map(int,archive['words']))
    families = {name:load_family(name,old,generator)
                for name in ('global','expanded66','expanded67_68')}
    assert [len(families[n]['new']) for n in families] == [1835,1803,1792]
    global_family = families['global']
    all_supports = old+global_family['new']
    search = best_four(all_supports)
    removed = sum(1 << c for c in FOUR)
    chosen = [s for s in all_supports if not(s & removed)]
    assert len(chosen) == 10261
    projection = [sum((row>>c&1)<<j for j,c in enumerate(FOUR)) for row in generator]
    assert len(rref(projection,4)[0]) == 4
    assert len(impose_zeros(generator,FOUR[2:])) == 26
    odd = odd_phases(global_family,removed)
    # Every support is in C^perp, so every dense restriction has even parity.
    assert all(bits(row&s)%2 == 0 for row in generator for s in chosen)
    assert geometry()['valid'] and not hex_geometry()[1]
    assert Q(14,32) <= Q(1,2) and quad_leq(0,Q(14,128),17)
    dimensions = {}
    for d,block_count in ((61,8555),(62,9145)):
        blocks = [b for b in sqs64() if max(b) < d]
        assert len(blocks) == block_count
        terms = dict(dense=(1<<26) if d==61 else 6*(1<<24),
                     weight16=len(chosen)*2048,weight4=block_count*16,poles=2*d)
        dimensions[str(d)] = dict(total=sum(terms.values()),terms=terms)
    # Validate the octad catalogue and prove completeness by MacWilliams' identity.
    with np.load(HERE.parent/'dim61_96/dual_sections.npz',allow_pickle=False) as archive:
        sample = [sum(1 << int(c) for c in row) for row in archive['section_supports']]
    g73 = sum(1<<e for e in EXPONENTS)
    parent = [g73<<i for i in range(36)]
    dual9 = sum(1<<c for c in DUAL_SUPPORT)
    assert all(bits(row&dual9)%2 == 0 for row in parent)
    head = [i for i in range(1,73) if i not in DUAL_SUPPORT]+[0]
    assert rref(keep_coordinates(impose_zeros(parent,DUAL_SUPPORT),head),64)[0] == rref(generator,64)[0]
    octads = json.loads((HERE/'dual_octads.json').read_text())
    macwilliams = sum(int(histogram[w])*sum((-1)**j*comb(w,j)*comb(64-w,8-j)
                      for j in range(9) if j<=w and 8-j<=64-w)
                      for w in range(65))
    assert macwilliams == 72*(1<<28)
    assert len(octads) == len(set(octads)) == 72 and all(bits(h)==8 for h in octads)
    assert all(bits(row&h)%2 == 0 for row in generator for h in octads)
    assert all(h in octads for h in sample)
    section_rows = []
    for h in octads:
        assert len(rref([local(h,row) for row in generator],8)[0]) == 7
        terms = dict(dense_section=70*(1<<21),
            original_weight16_section=sum(section_count(s,h,global_family['data']) for s in old),
            added_weight16_section=sum(section_count(s,h,global_family['data']) for s in global_family['new']),
            weight4_section=102900,ordinary_poles=112)
        section_rows.append((sum(terms.values()),h,terms))
    total,section,terms = max(section_rows,key=lambda row:(row[0],-row[1]))
    check_section_counts(global_family,old,section)
    # Check the independent affine four-support section count explicitly.
    four_section = 0
    for block in sqs64():
        m = sum(c < 8 for c in block)
        if m%2 == 0: four_section += comb(m,m//2)*(1<<(4-m))
    assert four_section == terms['weight4_section'] == 102900
    dimensions['63'] = dict(total=total,terms=terms)
    base = (1<<28)+30828*2048+10416*16+128
    terms64 = dict(dense=1<<28,weight16=(30828+1835)*2048,weight4=10416*16,poles=128)
    dimensions['64'] = dict(total=sum(terms64.values()),terms=terms64)
    # The independent BCH/affine-line owner proof is replayed by the 64--65D package.
    dimensions['65'] = dict(total=dimensions['64']['total']+68224+2,
        terms=dict(base64=dimensions['64']['total'],owner_gain=68224,new_poles=2))
    assert dimensions['64']['total'] == base+1835*2048
    folder = HERE.parent/'dim61_96'
    mother = np.fromfile(folder/'cyclic73.weight16.bin',dtype='<u8').reshape(-1,2)
    assert len(mother) == len(np.unique(mother,axis=0)) == 174324
    assert np.all(popcounts64(mother[:,0])+popcounts64(mother[:,1]) == 16)
    reduced,pivots = rref(parent,73)
    for lo,hi in mother:
        word = int(lo) | (int(hi)<<64)
        for row,pivot in zip(reduced,pivots):
            if word>>pivot&1: word ^= row
        assert not word
    heads = np.zeros(len(mother),dtype=np.uint64)
    for j,c in enumerate(head):
        heads |= ((mother[:,c//64]>>np.uint64(c%64))&np.uint64(1)) << np.uint64(j)
    choices = json.loads((folder/'shortening_support_search.json').read_text())['dimensions']
    netto = json.loads((HERE.parent/'continuation49_72/netto.json').read_text())
    expanded_checks = {}
    old_set = set(old)
    for d,name in ((66,'expanded66'),(67,'expanded67_68'),(68,'expanded67_68')):
        retained = np.ones(len(mother),dtype=bool)
        for c in choices[str(d)]['deleted_coordinates']:
            retained &= ((mother[:,c//64]>>np.uint64(c%64))&np.uint64(1)) == 0
        kept = heads[retained]
        ordinary = np.array([int(s) in old_set for s in kept])
        assert int(ordinary.sum()) == 30828
        extra = kept[~ordinary]
        assert len(kept) == choices[str(d)]['support_count']
        family = families[name]
        extra_coordinates = np.array([[(int(s)>>c)&1 for c in range(64)] for s in extra],dtype=np.float32)
        for a in family['new']:
            counts = popcounts64(extra & np.uint64(a))
            point = np.array([(a>>c)&1 for c in range(64)],dtype=np.float32)
            assert np.array_equal(counts,extra_coordinates @ point)
            assert int(counts.max()) <= 8
        expanded_checks[str(d)] = dict(family=name,extra_supports=len(extra),
            mixed_pairs_checked=len(extra)*len(family['new']))
        terms = dict(netto['dimensions'][str(d)]['terms'])
        terms['weight16'] += 2048*len(family['new'])
        dimensions[str(d)] = dict(total=sum(terms.values()),terms=terms)
    # These bounds cover every dense/sparse, sparse/four-support and axial pair.
    assert 2*16 == 2*4*4 == 32 and 2*8 <= 32
    for d,row in dimensions.items():
        row.update(previous_bound=PREVIOUS[int(d)],gain=row['total']-PREVIOUS[int(d)])
        assert row['total'] == sum(row['terms'].values()) and row['gain'] > 0
    report = dict(status='Verified: shared parity repairs and exact transfers to 61--68D',
        dimensions=dimensions,families={name:f['receipt'] for name,f in families.items()},
        parent64_words_enumerated=1<<28,
        parent64_weight_histogram={str(i):int(v) for i,v in enumerate(histogram) if v},
        four_zero_search=search,retained_weight16_supports=len(chosen),odd61=odd,
        section63=dict(selected_octad=section,candidates_checked=len(octads),
                       sign_counts_independently_enumerated=True),
        expanded_checks=expanded_checks,negative_control=negative_control(global_family),
        inherited_inputs='Existing [73,36,16] distance certificate and Netto designs; '
                        'the separate 64--65D verifier replays all 68,224 owners for 65D.',
        geometry='All shared bodies have norm64 and pair products <=32. '
                 '61D and 62D use the exact latitude and hexagonal maps with unit norm.')
    (HERE/'shared_layers.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__ == '__main__':
    main()
