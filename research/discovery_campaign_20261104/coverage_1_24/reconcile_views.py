"""Export active benchmark coverage views from the independently checked 24-row matrix.

The immutable incoming patch preserves the original author's audit, including its
historical unavailable-source observations; active views use the successful source.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from build_matrix import ROOT, HERE, check, sha, FIXTURES

BENCH = ROOT / 'discovery/benchmarks'
PATCH_SHA = '7ac4796cc14c476424e3a09888e499c077844828c9921b7508c8e81e384bd741'


def exports():
    matrix = check()
    rows = matrix['rows']
    provenance = json.loads((HERE / 'integration/provenance.json').read_text())
    if sha(ROOT / provenance['patch']) != PATCH_SHA:
        raise ValueError('Original incoming patch hash mismatch')
    if sha(BENCH / 'catalog.json') != provenance['catalogue_sha256_before_integration']:
        raise ValueError('Comparison catalogue changed during coverage integration')
    source = dict(rows[0]['source'])
    out = {'schema_version': 2, 'audit_date': '2026-11-04',
           'built_at_utc': matrix['built_at_utc'],
           'canonical_matrix': str((HERE / 'matrix.json').relative_to(ROOT)),
           'canonical_matrix_sha256': sha(HERE / 'matrix.json'),
           'incoming_patch': provenance,
           'scope': 'Current source-backed construction reproduction, not new records or worldwide novelty. Exactly 24 dimensions.',
           'status_semantics': {
               'imported': 'At least one exact target route has literal coordinate/Gram source ancestry.',
               'formula_generated': 'An exact target route constructs formulas at runtime, or a fixture has verified equality with the root formula generator.',
               'search_discovered': 'Historical autonomy is not certified by this coverage audit; separate held-out benchmarks record search results.',
               'search_execution': 'A target route executes a registered candidate search rather than only a fixed formula.',
               'exact_verified': 'At least one current target recipe passed independent exact replay.',
               'unsupported': 'No integrated exact target recipe in the audited inventory; generic operators and separate prototypes may exist.',
               'not_tested': 'No integrated target proof replay; this is not a failed search.'},
           'sources': {'cohn_successful_20261104': source}, 'rows': []}
    for row in rows:
        matched = [run for run in row['integrated_runs'] if run['cardinality'] == row['sourced_target']]
        fixtures = [run['recipe'] for run in matched if Path(run['recipe']).name in FIXTURES]
        out['rows'].append({
            'dimension': row['dimension'], 'cohn_lower_bound': row['sourced_target'],
            'bound_source': 'cohn_successful_20261104', 'source': row['source'],
            'target_kind': row['target_kind'], 'coverage_status': row['coverage_status'],
            'target_status': {'imported': any(run['origin'] in {'imported','mixed'} for run in matched),
                              'formula_generated': bool(fixtures) or any(run['origin'] in {'generated','mixed'} and not any(op.startswith('search.') for op in run['operators']) for run in matched),
                              'search_execution': any(any(op.startswith('search.') for op in run['operators']) for run in matched),
                              'search_discovered': False, 'exact_verified': bool(matched),
                              'unsupported': not bool(matched), 'not_tested': not bool(matched)},
            'recipe': row['implemented_recipe_or_none'],
            'verification_evidence': row['integrated_runs'],
            'formula_regeneration_verified_fixtures': fixtures,
            'numerical_search_recipes': row['implemented_numerical_search_recipes'],
            'separate_repository_evidence': row['separate_repository_evidence'],
            'remaining_gap': row['remaining_gap'], 'resource_limits': row['resource_limits'],
            'new_record': False})
    common = {'schema_version': 2, 'built_at_utc': matrix['built_at_utc'],
              'canonical_matrix': out['canonical_matrix'], 'canonical_matrix_sha256': out['canonical_matrix_sha256'],
              'source': source, 'incoming_patch_sha256': PATCH_SHA,
              'scope': 'Current local independent proof replay; original external timings and claims remain in preserved incoming patch; no new record.'}
    roots = dict(common, controls=[run for row in rows for run in row['integrated_runs'] if Path(run['recipe']).name in FIXTURES])
    other = dict(common, controls=[run for row in rows for run in row['integrated_runs'] if Path(run['recipe']).name not in FIXTURES])
    return {'coverage_1_24_20261104.json': out, 'root_controls_20261104.json': roots,
            'coverage_controls_20261104.json': other}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for name, value in exports().items():
        path = BENCH / name
        if args.check:
            if json.loads(path.read_text()) != value:
                raise ValueError(f'Stale active coverage view: {name}')
        else:
            path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'rows': 24, 'views': 3,
                      'mode': 'check' if args.check else 'write'}))


if __name__ == '__main__':
    main()
