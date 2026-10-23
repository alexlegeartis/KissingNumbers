#!/usr/bin/env python3
"""Read the construction ledger and recompute the five supplied current totals.

This driver does not verify angles or enumerate bodies; run verify.py for the
consistency replay. The missing latest 43D certificate cannot enter this list.
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
    assert {row['dimension'] for row in supplied} == {41, 42, 44, 45, 46}
    out = []
    for row in supplied:
        with np.load(ROOT / row['candidate'], allow_pickle=False) as model:
            body = sum(map(int, model['weights']))
            pure = sum(len(model[key]) for key in
                       ('pure_directions', 'pure_head_numerators') if key in model.files)
        assert body == row['body_points'] and pure == row['pure_points']
        assert body + pure == row['total']
        out.append(dict(row, population_scope=ledger()['population_scope']))
    return out


if __name__ == '__main__':
    for record in results():
        print(json.dumps(record))
