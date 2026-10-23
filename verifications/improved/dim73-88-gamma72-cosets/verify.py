"""Replay the nearestness certificates and all retained cross-head inequalities."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        replay(source, 'dim61_96/gamma72_coset_certificate.py')
        replay(source, 'dim61_96/gamma72_repaired_verify.py')
        for count in [4, 8, 16]:
            if (source / ('dim61_96/gamma72_scaled_%d.json' % count)).exists():
                replay(source, 'dim61_96/gamma72_scaled_verify_fast.py', '--seeds', str(count))
    print(json.dumps(dict(status='PASS', results=results())))


if __name__ == '__main__':
    main()
