#!/usr/bin/env python3
"""Read the construction ledger and recompute the six supplied current totals.

This driver does not verify angles or enumerate bodies; run verify.py for the
consistency replay, including the portable partial-family 43D certificate.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def ledger():
    return json.loads((ROOT / 'bounds.json').read_text(encoding='utf-8'))


def results():
    rows = ledger()['constructions']
    assert len({row['dimension'] for row in rows}) == len(rows)
    supplied = [row for row in rows if row['status'] == 'candidate_supplied']
    assert {row['dimension'] for row in supplied} == set(range(41, 47))
    out = []
    for row in supplied:
        with np.load(ROOT / row['candidate'], allow_pickle=False) as model:
            body = sum(map(int, model['weights']))
            pure = sum(len(model[key]) for key in
                       ('pure_directions', 'pure_head_numerators') if key in model.files)
        if 'explicit_points' in row:
            with np.load(ROOT / row['explicit_points'], allow_pickle=False) as points:
                explicit = len(points['body_indices'])
                assert explicit == row['explicit_body_points']
                body += explicit
        assert body == row['body_points'] and pure == row['pure_points']
        assert body + pure == row['total']
        out.append(dict(row, population_scope=ledger()['population_scope']))
    return out


if __name__ == '__main__':
    for record in results():
        print(json.dumps(record))
