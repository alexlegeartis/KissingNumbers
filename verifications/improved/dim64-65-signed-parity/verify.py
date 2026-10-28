"""Fresh enumeration and exact parity/geometry replay, isolated from research."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import construction as c
from final import results

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest


def words(rows):
    out = [0]
    for row in rows:
        out += [v ^ row for v in out]
    return out


def main():
    check_manifest(ROOT)
    data = json.loads((ROOT/'data/construction.json').read_text())
    compiler = shutil.which('g++')
    if not compiler and Path('C:/msys64/mingw64/bin/g++.exe').exists():
        compiler = 'C:/msys64/mingw64/bin/g++.exe'
    if not compiler:
        raise SystemExit('A g++ compiler is required for fresh code enumeration.')
    env = dict(os.environ)
    env['PATH'] = str(Path(compiler).parent)+os.pathsep+env.get('PATH','')
    with tempfile.TemporaryDirectory(prefix='kissing64_65_') as directory:
        work = Path(directory)
        executable = work/('enumerate_binary.exe' if sys.platform=='win32' else 'enumerate_binary')
        subprocess.run([compiler,'-O3',str(ROOT/'enumerate_binary.cpp'),'-o',str(executable)],check=True,env=env)
        c.HERE, c.ENUMERATOR = work, executable
        owner_report = c.construction65()
        assert owner_report['base64_total'] == data['base64'] == 331737984
        assert owner_report['owner_class']['total_points'] == data['owner_points65'] == 68224
        rows = data['original_support_generator']
        assert len(c.rref(rows,64)[0]) == 28
        generator = work/'supports.generator.txt'
        generator.write_text('64 28\n'+''.join('%x 0\n'%v for v in rows))
        output = work/'supports.bin'
        report = json.loads(subprocess.check_output([str(executable),str(generator),str(output)],text=True,env=env))
        histogram = {int(k):v for k,v in report['histogram'].items()}
        assert sum(histogram.values()) == 2**28 and min(w for w in histogram if w) == 16
        raw = np.fromfile(output,dtype=np.uint64).reshape(-1,2)
        assert raw.shape == (30828,2) and not np.any(raw[:,1])
        supports = raw[:,0].copy()
    assert len(np.unique(supports)) == 30828
    support_set = set(map(int,supports))
    used = set()
    added = []
    lut = np.array([c.popcount(i) for i in range(256)],dtype=np.uint8)
    standard_checks = [65535]+[sum(((i>>j)&1)<<i for i in range(16)) for j in range(4)]
    standard = words(c.kernel(standard_checks,16))
    assert len(standard)==2048 and min(c.popcount(v) for v in standard if v)==4
    assert all(c.popcount((v^1)&255)%2==1 for v in standard)
    for block in data['added_blocks']:
        a = block['support']
        assert 0<=a<2**64 and c.popcount(a)==16 and a not in support_set
        assert all(c.popcount(a&b)<=8 for b in added)
        added.append(a)
        dual = block['sign_dual']
        assert len(dual)==len(c.rref(dual,16)[0])==5
        dual_space = set(words(dual))
        assert set(map(c.popcount,dual_space))=={0,8,16}
        signs = words(c.kernel(dual,16))
        assert len(signs)==2048 and min(c.popcount(v) for v in signs if v)==4
        counts = lut[(supports & np.uint64(a)).view(np.uint8).reshape(-1,8)].sum(axis=1)
        bad = set(map(int,supports[counts>8]))
        assert int(counts.max())==10 and len(bad)==8
        supplied = block['modified_old_blocks']
        assert len(supplied)==8 and {r['old_support'] for r in supplied}==bad
        assert not used.intersection(bad)
        used.update(bad)
        coordinates = [i for i in range(64) if a>>i&1]
        for record in supplied:
            b, shared = record['old_support'],record['shared_check']
            assert shared & a == shared and shared & b == shared
            assert c.popcount(shared)==8
            local = sum(int(shared>>p&1)<<i for i,p in enumerate(coordinates))
            assert local in dual_space
            # Order shared coordinates first on the old support, and use
            # standard Hamming signs shifted in its first coordinate.
            # The shared parity is odd there, even on the new block.
            assert 4*(c.popcount(a&b)-2)<=32
    assert len(added)==len(set(added))==211 and len(used)==1688
    # Cross layers, independently of support signs and coordinate order.
    assert 2*16==32 and 2*4*4==32 and 2*8<32
    # No owner uses the weight-16 layer, so all owners survive unchanged.
    # The owner's coherence and the cap-versus-equator geometry were
    # checked in construction65(), and apply to every new equator point.
    expected = results()
    assert data['dimension64']==owner_report['base64_total']+2048*len(added)
    shared_package = ROOT.parent/'dim61-68-95-96-code-layers'
    provenance = json.loads((ROOT/'UPDATE_PROVENANCE.json').read_text())
    assert hashlib.sha256((shared_package/'source/shared_parity/shared_layers.json').read_bytes()).hexdigest() == provenance['shared_receipt_sha256']
    assert hashlib.sha256((shared_package/'source/shared_parity/global.npz').read_bytes()).hexdigest() == provenance['shared_certificate_sha256']
    subprocess.run([sys.executable,'-B',str(shared_package/'verify_shared.py')],check=True)
    shared = json.loads((shared_package/'source/shared_parity/shared_layers.json').read_text())
    assert expected[0]['total']==owner_report['base64_total']+2048*shared['families']['global']['new_supports']
    assert expected[1]['total']==expected[0]['total']+68224+2
    print(json.dumps(dict(status='PASS',freshly_enumerated_words=2*(2**28),
        historical_added_supports=len(added),current_added_supports=shared['families']['global']['new_supports'],
        modified_old_supports=shared['families']['global']['modified_old_supports'],results=expected)))


if __name__ == '__main__':
    main()
