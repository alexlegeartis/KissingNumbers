"""Expose the certified bound to the repository audit."""
import json
from pathlib import Path


def results():
    here = Path(__file__).resolve().parent
    w = json.loads((here / 'witness.json').read_text())
    r = json.loads((here / 'verification.recorded.json').read_text())
    assert len(r['sizes']) == len(w['families']) == r['families'] == 53
    assert r['D4_calibrations'] == w['calibrations'] and r['cross_family_tests'] == 1378
    assert sum(r['sizes']) + len(w['pure_head_signs']) == r['points'] == w['total']
    return [{'dimension': w['dimension'], 'total': w['total']}]


if __name__ == '__main__':
    print(json.dumps(results()))
