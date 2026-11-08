"""Full read-only finite replay of the reconstructed 49–60D P48 certificates."""
from pathlib import Path
import hashlib
import json
import sys
from threadpoolctl import threadpool_limits
import final
import independent_audit
import tail_and_lift
import verify_heads

D = Path(__file__).resolve().parent


def verify_manifest(root=D):
    root = Path(root).resolve()
    manifest = json.loads((root / 'MANIFEST.json').read_text())
    for name, item in manifest['files'].items():
        path = (root / name).resolve()
        if root not in path.parents:
            raise ValueError('Manifest path leaves package')
        raw = path.read_bytes()
        if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise ValueError('Manifest mismatch: ' + name)
    return len(manifest['files'])


def main():
    if not __debug__:
        raise SystemExit('Run without -O: the independent checkers require assertions.')
    manifest_digest = hashlib.sha256((D / 'MANIFEST.json').read_bytes()).hexdigest()
    manifest_files = verify_manifest()
    metadata = json.loads((D / 'certificate.json').read_text())
    derived = final.results()
    tails = tail_and_lift.verify(D / 'source')
    base, gram, base_pairs = verify_heads.source(D / 'source')
    reports = []
    for row, count, tail in zip(metadata['rows'], derived, tails['tails']):
        path = D / ('heads%d.npz' % row['dimension'])
        local = verify_heads.certificate(path, base, gram, row['classes'])
        global_check = independent_audit.main(path, row['poles'], D / 'source')
        if (tail['ambient_dimension'] != row['dimension']
                or tail['antipodal_pairs'] != row['classes']
                or tail['pole_points'] != row['poles']
                or local['selected_lines'] != row['selected_lines']
                or global_check['selected_lines'] != row['selected_lines']
                or local['certificate_sha256'] != global_check['certificate_sha256']
                or global_check['conditional_kissing_count'] != count['total']):
            raise ValueError('Independent head, tail and count replays disagree')
        report = dict(dimension=row['dimension'], lower_bound=count['total'],
                      gain=count['total'] - row['previous_repository_bound'],
                      classes=row['classes'], selected_lines=local['selected_lines'],
                      extra_pairs=local['extra_pairs'], certificate_sha256=local['certificate_sha256'],
                      int64_local_frame='PASS', exact_binary64_global_frame='PASS')
        reports.append(report)
        print(json.dumps(report), flush=True)
    if (verify_manifest() != manifest_files or
            hashlib.sha256((D / 'MANIFEST.json').read_bytes()).hexdigest() != manifest_digest):
        raise ValueError('Manifest changed during replay')
    result = dict(status='PASS', rows=reports, base_pairs=base_pairs,
                  tail_pairs=sum(t['pairs_checked'] for t in tails['tails']),
                  independent_pole_pairs=tails['independent_poles']['pairs_checked'],
                  manifest_files=manifest_files, lift_pair_categories=tails['lift'],
                  inherited_inputs='Published identification of P48 Gram lattice, minimum squared norm six and 52416000 minimal vectors; no new shell enumeration.',
                  scope='Two exact finite head proofs, complete direction/pole pair checks and algebraic lift proof. Not an expanded all-pairs test of the full lifted codes; no global-optimality or worldwide-priority claim.')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    with threadpool_limits(limits=2):
        main()
