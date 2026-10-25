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
    for dim, row in mixed['dimensions'].items():
        assert row['bound'] == sum(row['terms'].values())
        for registered in rows:
            if registered['dimension'] == int(dim):
                registered['total'] = max(registered['total'], row['bound'])
    return sorted(rows, key=lambda row: row['dimension'])


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
