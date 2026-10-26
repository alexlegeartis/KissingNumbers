"""Two exact replays of the current 95--96D support certificate."""
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results

def validate_report(source):
    fresh = json.loads((source/'dim61_96/mixed_pair_certificate.json').read_text())
    recorded = json.loads((ROOT/'source/dim61_96/mixed_pair_certificate.json').read_text())
    assert {k:v for k,v in fresh.items() if k != 'seconds'} == {
        k:v for k,v in recorded.items() if k != 'seconds'}
    bits = json.loads((source/'dim61_96/support_bits_verification.json').read_text())
    assert bits == json.loads((ROOT/'source/dim61_96/support_bits_verification.json').read_text())
    for key in bits.keys() - {'status'}:
        assert bits[key] == fresh[key], key
    expected = {r['dimension']:r['total'] for r in results() if r['dimension'] in (95, 96)}
    assert {int(d):r['bound'] for d,r in fresh['dimensions'].items()} == expected
    return fresh


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        replay(source, 'dim61_96/mixed_pair_certificate.py')
        replay(source, 'dim61_96/verify_support_bits.py')
        report = validate_report(source)
    print(json.dumps({'status':'PASS', 'dimensions':report['dimensions'],
                      'support_pairs_checked_twice':report['support_pairs_checked']}))


if __name__ == '__main__':
    main()
