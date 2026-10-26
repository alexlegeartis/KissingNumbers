"""Derive the current code-layer totals from the supplied certificates."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / 'source' / name).read_text(encoding='utf-8'))


def results():
    rows = []
    counts = read('dim49plus/latitude_compression.json')['counts']
    assert counts['total'] == counts['dense'] + counts['weight15'] + counts['weight2']
    rows.append(dict(dimension=61, total=counts['total']))
    for name in ['dual_sections.json', 'expanded_support_layers.json']:
        report = read('dim61_96/' + name)
        for dim, row in report['dimensions'].items():
            assert row['total'] == sum(row['terms'].values())
            rows.append(dict(dimension=int(dim), total=row['total']))
    section = read('dim61_96/section95_search.json')
    best = section[section['best_method']]
    assert sum(best['terms'].values()) == section['dimension95_total'] == best['total']
    rows.append(dict(dimension=95, total=best['total']))
    row = read('dim61_96/support_doubling.json')['dimension96']
    assert row['total'] == sum(row['terms'].values())
    rows.append(dict(dimension=96, total=row['total']))
    pairing = read('dim61_96/support_pairing_certificate.json')
    assert pairing['status'].startswith('Verified:')
    for dim, row in pairing['dimensions'].items():
        assert row['bound'] == sum(row['terms'].values())
        for registered in rows:
            if registered['dimension'] == int(dim):
                registered['total'] = max(registered['total'], row['bound'])
    mixed = read('dim61_96/mixed_pair_certificate.json')
    assert mixed['status'].startswith('Verified:')
    assert mixed['support_pairs_checked'] == mixed['supports']*(mixed['supports']-1)//2
    assert sum(mixed['section_intersection_histogram']) == mixed['supports']
    assert mixed['dimensions']['96']['terms']['new_support24'] == mixed['supports']*16384
    assert mixed['dimensions']['95']['terms']['new_support24_section'] == sum(
        count*signs for count, signs in zip(mixed['section_intersection_histogram'],
                                          mixed['retained_sign_counts']))
    for dim, row in mixed['dimensions'].items():
        assert row['bound'] == sum(row['terms'].values())
        assert row['gain'] == row['bound']-row['previous_bound']
        for registered in rows:
            if registered['dimension'] == int(dim):
                registered['total'] = max(registered['total'], row['bound'])
    return sorted(rows, key=lambda row: row['dimension'])


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
