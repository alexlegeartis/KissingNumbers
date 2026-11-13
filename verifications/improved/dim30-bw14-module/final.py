"""Read the registered count from the complete coordinate certificate."""
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parent


def results():
    with np.load(ROOT / 'data/coordinates30_248020.npz', allow_pickle=False) as z:
        n, d = z['points'].shape
    if (n, d) != (248020, 30):
        raise ValueError('Certificate dimensions')
    return [dict(dimension=d, total=n)]


if __name__ == '__main__':
    for row in results():
        print(json.dumps(row))
