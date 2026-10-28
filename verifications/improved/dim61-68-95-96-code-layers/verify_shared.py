"""Replay the shared parity certificates; optionally re-enumerate all 2^36 words."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results


def validate(source):
    relative = Path('shared_parity/shared_layers.json')
    fresh = json.loads((source/relative).read_text())
    assert fresh == json.loads((ROOT/'source'/relative).read_text())
    current = {row['dimension']:row['total'] for row in results()}
    totals = {int(d):row['total'] for d,row in fresh['dimensions'].items() if int(d) in current}
    assert totals == {d:current[d] for d in totals}
    return totals


def full_codes(source):
    compiler = shutil.which('g++')
    if not compiler: raise SystemExit('--full-codes requires g++.')
    folder = source/'dim61_96'
    executable = folder/('shared_enumerator.exe' if sys.platform=='win32' else 'shared_enumerator')
    subprocess.run([compiler,'-O3','-mpopcnt','-mbmi',str(folder/'enumerate_binary.cpp'),
                    '-o',str(executable)],check=True)
    sys.path.insert(0,str(source/'dim49plus'))
    from binary_62_code import EXPONENTS
    g = sum(1<<e for e in EXPONENTS)
    generator = folder/'rebuilt73.generator.txt'
    generator.write_text('73 36\n'+''.join(f'{(g<<i)&((1<<64)-1):x} {(g<<i)>>64:x}\n' for i in range(36)))
    actual = json.loads(subprocess.check_output([str(executable),str(generator)],text=True))
    recorded = json.loads((folder/'cyclic73.json').read_text())
    assert (actual['n'],actual['k']) == (73,36)
    assert actual['histogram'] == recorded['histogram']
    print('PASS: independent C++ enumeration checked all 68,719,476,736 mother-code words.',flush=True)


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        # Recheck every inherited compression and cross-layer inequality too.
        replay(source,'significant49_72/compressed_layers.py')
        replay(source,'shared_parity/shared_layers.py')
        totals = validate(source)
        if '--full-codes' in sys.argv: full_codes(source)
    print(json.dumps(dict(status='PASS',dimensions=totals,full_codes='--full-codes' in sys.argv)))


if __name__ == '__main__':
    main()
