"""Formula-generated root-system controls through the existing integer operator.

This is benchmark recipe generation, not a new production operator or search.
No record coordinates are read. See docs/COVERAGE_1_24.md for proof and sources.
"""
import argparse
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import platform
import time

from discovery.core import RunContext, digest
from discovery.runner import output_artifact, run_program


ROOT_NAMES = {1: 'a1', 2: 'a2', 3: 'd3', 4: 'd4', 5: 'd5',
              6: 'e6', 7: 'e7', 8: 'e8'}
EXPECTED_COUNTS = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240}


def root_recipe(dimension):
    """Return a deterministic v1 recipe with exact integer coordinates.

    E8 is scaled by two: the D8 part has two entries +/-2 and the other
    part has eight entries +/-1 with an even number of minus signs.
    E7 and E6 are its intersections with x6+x7=0 and additionally x5+x7=0
    (zero-based indices). Their ranks are checked, never assumed.
    """
    if type(dimension) is not int or dimension not in ROOT_NAMES:
        raise ValueError('Root benchmark dimension must be an integer from 1 to 8')
    if dimension == 1:
        rows, norm = [[-1], [1]], 1
    elif dimension == 2:
        rows = []
        for i, j in combinations(range(3), 2):
            for sign in (-1, 1):
                row = [0] * 3
                row[i], row[j] = sign, -sign
                rows.append(row)
        norm = 2
    else:
        ambient = dimension if dimension <= 5 else 8
        scale = 1 if dimension <= 5 else 2
        rows = []
        for i, j in combinations(range(ambient), 2):
            for a, b in product((-scale, scale), repeat=2):
                row = [0] * ambient
                row[i], row[j] = a, b
                rows.append(row)
        if dimension >= 6:
            rows.extend(list(signs) for signs in product((-1, 1), repeat=8)
                        if sum(x < 0 for x in signs) % 2 == 0)
            if dimension <= 7:
                rows = [x for x in rows if x[6] + x[7] == 0]
            if dimension == 6:
                rows = [x for x in rows if x[5] + x[7] == 0]
        norm = 2 * scale**2
    name = ROOT_NAMES[dimension]
    return {
        'schema_version': 1,
        'name': f'roots_{name}',
        'description': 'Known root-system control generated from formulas; no new SOTA claim.',
        'provenance': {
            'method': 'formula-generated benchmark coordinates; no search or coordinate import',
            'generator': 'discovery/benchmarks/root_controls.py',
            'construction_reference': 'https://arxiv.org/html/2411.04916v3#S1',
            'root_system': name.upper(),
        },
        'nodes': [{'id': 'roots', 'op': 'geometry.integer_code',
                   'params': {'coordinates': sorted(rows), 'norm_squared': norm}}],
        'outputs': ['roots'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--dimensions', type=int, nargs='+', choices=range(1, 9),
                        default=list(range(1, 9)))
    args = parser.parse_args()
    if len(set(args.dimensions)) != len(args.dimensions):
        parser.error('Dimensions must be distinct')
    if args.outdir.exists():
        parser.error('Use a new output directory; benchmark files are not overwritten')
    args.outdir.mkdir(parents=True)
    context = RunContext(args.outdir)
    controls = []
    for dimension in args.dimensions:
        recipe = root_recipe(dimension)
        name = recipe['name']
        start = time.monotonic()
        receipt = run_program(recipe, context)
        artifact = output_artifact(receipt)
        if (artifact.verification.status != 'exact'
                or artifact.facts['dimension'] != dimension
                or artifact.facts['cardinality'] != EXPECTED_COUNTS[dimension]):
            raise RuntimeError(f'{name} failed its independently verified benchmark expectation')
        for suffix, value in (('recipe', recipe), ('receipt', receipt)):
            (args.outdir / f'{name}.{suffix}.json').write_text(json.dumps(value, indent=2) + '\n')
        row = {'name': name, 'generation': 'formula-generated',
               'verification': 'exact-verified', **artifact.facts,
               'program_sha256': digest(recipe),
               'seconds': round(time.monotonic() - start, 6),
               'recipe': f'{name}.recipe.json', 'receipt': f'{name}.receipt.json',
               'independent_terminal_replay': True, 'new_sota': False}
        controls.append(row)
        print(json.dumps(row), flush=True)
    summary = {
        'schema_version': 1, 'python': platform.python_version(),
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source_sha256': receipt['source_sha256'],
        'budget': context.budget, 'controls': controls,
        'scope': 'Known formula-generated root controls; measured on this executor; no search or novelty claim.',
    }
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
