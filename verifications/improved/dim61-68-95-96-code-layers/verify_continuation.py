"""Complete exact replays of the new 61--63D and 66--68D constructions."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results


def validate(source):
    totals = {}
    for name in ['bridge61', 'signed_sections', 'netto']:
        fresh = json.loads((source/'continuation49_72'/(name+'.json')).read_text())
        recorded = json.loads((ROOT/'source/continuation49_72'/(name+'.json')).read_text())
        assert fresh == recorded, name
        if name == 'bridge61':
            totals[61] = fresh['bound']
        else:
            totals.update({int(d): row['total'] for d, row in fresh['dimensions'].items()})
    expected = {row['dimension']:row['total'] for row in results() if row['dimension'] in totals}
    assert totals == expected
    return totals


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        for name in ['bridge61', 'signed_sections', 'netto']:
            replay(source, 'continuation49_72/'+name+'.py')
        totals = validate(source)
    print(json.dumps(dict(status='PASS', dimensions=totals)))


if __name__ == '__main__':
    main()
