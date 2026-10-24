"""Check package integrity and replay the exact 28D verifier in isolation."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, temporary_source, replay
from final import results


if __name__ == '__main__':
    check_manifest(ROOT)
    print(results(), flush=True)
    with temporary_source(ROOT) as source:
        replay(source, 'verify.py', 'D28-205232.txt')
