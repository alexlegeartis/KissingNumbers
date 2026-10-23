"""Replay the binary compression, code layers and central sections in isolation."""
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        replay(source, 'dim49plus/latitude_compression.py')
        replay(source, 'dim61_96/weighted_layers.py')
        replay(source, 'dim61_96/dual_sections.py')
        if '--full-codes' in sys.argv:
            compiler = shutil.which('g++')
            if not compiler and Path('C:/msys64/mingw64/bin/g++.exe').exists():
                compiler = 'C:/msys64/mingw64/bin/g++.exe'
            if not compiler:
                raise SystemExit('--full-codes requires a C++ compiler (g++).')
            folder = source / 'dim61_96'
            executable = folder / ('enumerate_binary.exe' if sys.platform=='win32' else 'enumerate_binary')
            subprocess.run([compiler, '-O3', '-mpopcnt', '-mbmi', str(folder/'enumerate_binary.cpp'),
                            '-o', str(executable)], check=True)
            output = subprocess.check_output([str(executable), str(folder/'cyclic73.generator.txt')], text=True)
            actual = json.loads(output)
            expected = json.loads((folder/'cyclic73.json').read_text())
            assert (actual['n'], actual['k']) == (73, 36)
            assert actual['histogram'] == expected['histogram']
            print('PASS: all 68,719,476,736 binary codewords checked.', flush=True)
        replay(source, 'dim61_96/expanded_support_layers.py')
        replay(source, 'dim61_96/support_doubling.py')
        replay(source, 'dim61_96/section95_search.py')
        fresh = {}
        report = json.loads((source/'dim49plus/latitude_compression.json').read_text())
        fresh[61] = report['counts']['total']
        for name in ['dual_sections.json','expanded_support_layers.json']:
            report = json.loads((source/'dim61_96'/name).read_text())
            fresh.update({int(d):r['total'] for d,r in report['dimensions'].items()})
        fresh[95] = json.loads((source/'dim61_96/section95_search.json').read_text())['dimension95_total']
        fresh[96] = json.loads((source/'dim61_96/support_doubling.json').read_text())['dimension96']['total']
        assert fresh == {r['dimension']:r['total'] for r in results()}
    print(json.dumps(dict(status='PASS', results=results(), full_codes='--full-codes' in sys.argv)))


if __name__ == '__main__':
    main()
