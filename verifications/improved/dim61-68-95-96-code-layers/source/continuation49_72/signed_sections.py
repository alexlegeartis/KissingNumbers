"""Exact 62--63D central sections of the parity-augmented 64D construction.

All sign orders and coset shifts are explicit. Packed counts are compared
with expanded integer signs. No optimization score is used as proof.
"""
import hashlib
import json
from math import comb
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'dim61_96'))
from weighted_layers import kernel, packed, rref, sqs64
from dual_sections import dense_section_count


def binary_words(basis):
    result = [0]
    for row in basis:
        result += [word ^ row for word in result]
    return result


def count_bits(words, masks):
    sizes = [mask.bit_count() for mask in masks]
    if any(size % 2 for size in sizes):
        return 0
    return sum(all((word & mask).bit_count() == size//2
                   for mask, size in zip(masks, sizes)) for word in words)


def count_coordinates(words, masks):
    signs = np.array([[1-2*((word >> i) & 1) for i in range(16)]
                      for word in words], dtype=np.int8)
    retained = np.ones(len(words), dtype=bool)
    for mask in masks:
        columns = [i for i in range(16) if mask >> i & 1]
        retained &= signs[:, columns].sum(axis=1, dtype=np.int16) == 0
    return int(retained.sum())


def checked_count(words, masks):
    value = count_bits(words, masks)
    assert value == count_coordinates(words, masks)
    return value


def ordered_masks(coordinates, sections):
    return [sum(int(section >> c & 1) << i for i, c in enumerate(coordinates))
            for section in sections]


def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are required.')
    parent_path = HERE / 'parent64.json'
    parent = json.loads(parent_path.read_text())
    shifts = json.loads((HERE / 'section_shifts.json').read_text())
    folder = HERE.parent / 'dim61_96'
    with np.load(folder / 'code64_supports.npz', allow_pickle=False) as archive:
        generator, supports = archive['generator'].copy(), list(map(int, archive['words']))
    assert len(supports) == len(set(supports)) == 30828
    assert all(word.bit_count() == 16 for word in supports)
    assert rref(parent['original_support_generator'], 64)[0] == rref(packed(generator), 64)[0]
    rows, pivots = rref(packed(generator), 64)
    assert len(rows) == 28
    for word in supports:
        for row, pivot in zip(rows, pivots):
            if word >> pivot & 1:
                word ^= row
        assert word == 0
    with np.load(folder / 'dual_sections.npz', allow_pickle=False) as archive:
        section_coordinates = archive['section_supports'].copy()
    all_sections = [sum(1 << int(c) for c in row) for row in section_coordinates]
    assert all_sections[0] & all_sections[1] == 0
    assert all(section.bit_count() == 8 for section in all_sections)
    assert all((row & section).bit_count() % 2 == 0
               for row in packed(generator) for section in all_sections)

    checks = [65535]+[sum(((i >> j) & 1) << i for i in range(16)) for j in range(4)]
    standard = binary_words(kernel(checks, 16))
    assert len(standard) == len(set(standard)) == 2048 and min(word.bit_count() for word in standard if word) == 4
    assert all((word & 255).bit_count() % 2 == 0 for word in standard)
    additions, used = [], set()
    for block in parent['added_blocks']:
        support = block['support']
        assert support.bit_count() == 16 and support not in supports
        assert all((support & other).bit_count() <= 8 for other in additions)
        additions.append(support)
        bad = {old for old in supports if (old & support).bit_count() > 8}
        assert len(bad) == 8 and all((old & support).bit_count() == 10 for old in bad)
        assert bad == {record['old_support'] for record in block['modified_old_blocks']}
        assert not used.intersection(bad)
        used.update(bad)
        sign_words = binary_words(kernel(block['sign_dual'], 16))
        assert len(sign_words) == len(set(sign_words)) == 2048 and min(word.bit_count() for word in sign_words if word) == 4
    assert len(additions) == 211 and len(used) == 1688

    old_report = json.loads((folder / 'dual_sections.json').read_text())
    dimensions = {}
    for dimension in (62, 63):
        sections = all_sections[:64-dimension]
        assert len(shifts[str(dimension)]) == len(parent['added_blocks'])
        cache = {}

        def ordinary_count(support):
            shape = tuple(sorted(((support & section).bit_count() for section in sections), reverse=True))
            if shape not in cache:
                offset, masks = 0, []
                for size in shape:
                    masks.append(((1 << size)-1) << offset)
                    offset += size
                cache[shape] = checked_count(standard, masks)
            return cache[shape]

        original_sparse = sum(ordinary_count(support) for support in supports)
        assert original_sparse == old_report['dimensions'][str(dimension)]['terms']['weight16_section']
        old_before = old_after = new_points = 0
        for block, shift in zip(parent['added_blocks'], shifts[str(dimension)]):
            assert 0 <= shift < 2**16
            support = block['support']
            coordinates = [c for c in range(64) if support >> c & 1]
            new_words = [word ^ shift for word in binary_words(kernel(block['sign_dual'], 16))]
            new_points += checked_count(new_words, ordered_masks(coordinates, sections))
            for record in block['modified_old_blocks']:
                old, shared = record['old_support'], record['shared_check']
                assert shared & support == shared and shared & old == shared and shared.bit_count() == 8
                local = sum(int(shared >> c & 1) << i for i, c in enumerate(coordinates))
                new_parity = (shift & local).bit_count() % 2
                assert all((word & local).bit_count() % 2 == new_parity for word in new_words)
                # Put the shared eight-set first, then group by the disjoint sections.
                order = []
                for inside in (True, False):
                    for label in range(len(sections)+1):
                        order += [c for c in range(64) if old >> c & 1
                                  and bool(shared >> c & 1) == inside
                                  and next((i for i, section in enumerate(sections) if section >> c & 1), len(sections)) == label]
                assert len(order) == len(set(order)) == 16
                assert sum(1 << c for c in order[:8]) == shared
                target_parity = 1-new_parity
                old_words = [word ^ target_parity for word in standard]
                assert all((word & 255).bit_count() % 2 != new_parity for word in old_words)
                # Opposite parity on the shared subset forces a sign disagreement.
                assert 4*((support & old).bit_count()-2) <= 32
                old_before += ordinary_count(old)
                old_after += checked_count(old_words, ordered_masks(order, sections))

        dense, rank = dense_section_count(generator,
            [[c for c in range(64) if section >> c & 1] for section in sections])
        assert dense == old_report['dimensions'][str(dimension)]['terms']['dense_section']
        # Relabel the independent affine four-support layer so the section sets
        # are consecutive affine 3-spaces. All cross-layer bounds ignore labels.
        sparse4 = 0
        for block in sqs64():
            shape = [sum(8*i <= c < 8*(i+1) for c in block) for i in range(len(sections))]
            if all(size % 2 == 0 for size in shape):
                count = 2**(4-sum(shape))
                for size in shape:
                    count *= comb(size, size//2)
                sparse4 += count
        poles = 2*(64-8*len(sections))
        terms = dict(dense_section=dense, original_weight16_section=original_sparse-old_before+old_after,
                     added_weight16_section=new_points, weight4_section=sparse4, ordinary_poles=poles)
        previous = old_report['dimensions'][str(dimension)]['total']
        total = sum(terms.values())
        assert total > previous
        dimensions[str(dimension)] = dict(total=total, previous_bound=previous, gain=total-previous,
            terms=terms, section_supports=[[c for c in range(64) if s >> c & 1] for s in sections],
            projection_rank=rank, repaired_old_before=old_before, repaired_old_after=old_after)
    report = dict(status='Verified: explicit parity-compatible parent and two exact sign-count representations',
                  dimensions=dimensions, parent_additions=211, disjoint_repaired_old_supports=1688,
                  parent_certificate_sha256=hashlib.sha256(parent_path.read_bytes()).hexdigest(),
                  source_scope='The original [64,28,16] code and distance proof are unchanged. All new supports, parity repairs, coset shifts and section counts are checked here.',
                  geometry='Every section is central in a norm64 kissing configuration; orthogonal identification preserves norms and inner products.')
    (HERE / 'signed_sections.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
