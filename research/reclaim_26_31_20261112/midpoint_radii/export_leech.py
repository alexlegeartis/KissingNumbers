"""Regenerate the complete norm-32 Leech shell from the repository generator."""
from pathlib import Path
import argparse
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'source_and_layers'))
from source_paths import fresh_leech_shell

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    points = fresh_leech_shell()
    if points.shape != (196560, 24):
        raise ValueError('Leech shell population')
    if not np.all(np.sum(points.astype(np.int64)**2, axis=1) == 32):
        raise ValueError('Leech shell norms')
    np.save(args.output, points)
    print(f'Wrote {len(points)} norm-32 vectors to {args.output}')
