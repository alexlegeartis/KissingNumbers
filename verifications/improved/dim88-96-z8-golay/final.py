"""Read the verified Z8-lift sections; all counts have two exact replays."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def results():
    certificate = json.loads((ROOT/'data/replay.json').read_text())
    recorded = json.loads((ROOT/'bounds.json').read_text())
    assert certificate['status'] == recorded['status'] == 'PASS'
    assert certificate['z8']['binary_union_words'] == 2**37
    assert certificate['z8']['low_lee_coset_exclusions'] == 759*5489
    assert certificate['support']['support_pairs_checked'] == 26022*26021//2
    rows = []
    for row in recorded['rows']:
        dimension = row['dimension']
        source = certificate['dimensions'][str(dimension)]
        assert source['total'] == row['total'] == sum(row['terms'].values())
        assert source['terms'] == row['terms']
        assert len(source['constraint_rows']) == 96-dimension
        assert row['gain'] == row['total']-row['previous_repository'] > 0
        rows.append(dict(dimension=dimension, total=row['total']))
    assert [r['dimension'] for r in rows] == list(range(88,97))
    return rows


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
