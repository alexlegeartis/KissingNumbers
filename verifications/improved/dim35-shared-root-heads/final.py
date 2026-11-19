"""Derive the registered size from the compact construction."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent

def results():
    source = json.loads((ROOT / 'code_data.json').read_text())
    with np.load(ROOT / 'shared_positive.npz', allow_pickle=False) as z:
        stars, heads, cut = (len(z[k]) for k in ('new_stars', 'heads', 'removed_roots'))
    total = (1 << len(source['dense_generator'])) + len(source['full_parity_supports']) * 128
    total += 2 * 35 * 34 - cut + stars + heads
    return [{'dimension': 35, 'total': total}]

if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
