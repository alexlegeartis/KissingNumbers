"""Reproduce the finite support search; verify the resulting NPZ separately.

The supplied starting checkpoint already includes two insertions from an
earlier walk. The default run continues it using deterministic integer tests.
It does not train or reproduce PackingStar's reinforcement-learning system.
"""
import argparse
import hashlib
from pathlib import Path
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE / 'start_supports.npz')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=282029)
    parser.add_argument('--trials', type=int, default=100000000)
    parser.add_argument('--neutral-period', type=int, default=256)
    parser.add_argument('--compiler', default='g++')
    args = parser.parse_args()
    if args.trials < 0 or args.neutral_period < 1 or args.seed < 0:
        parser.error('trials/seed must be nonnegative and neutral-period positive')
    if args.output.resolve() == args.input.resolve():
        parser.error('output must differ from the input checkpoint')
    with np.load(args.input, allow_pickle=False) as data:
        low, high, sign = (data[k] for k in
                           ('supports_low', 'supports_high', 'sign_generator'))
    with tempfile.TemporaryDirectory(prefix='kissing_fill_correct_') as directory:
        work = Path(directory)
        source, result, executable = work/'input.txt', work/'output.txt', work/'search'
        with source.open('w') as stream:
            stream.write(str(len(low))+'\n')
            stream.writelines(f'{int(a)} {int(b)}\n' for a,b in zip(low,high))
        subprocess.run([args.compiler, '-O3', '-mpopcnt', '-std=c++17',
                        str(HERE/'fill_correct.cpp'), '-o', str(executable)], check=True)
        subprocess.run([str(executable), str(source), str(result), str(args.seed),
                        str(args.trials), str(args.neutral_period)], check=True)
        rows = np.loadtxt(result, skiprows=1, dtype=np.uint64, ndmin=2)
        assert rows.shape == (int(result.read_text().splitlines()[0]), 2)
        np.savez_compressed(args.output, supports_low=rows[:,0],
                            supports_high=rows[:,1], sign_generator=sign)
    print('Candidate:', args.output)
    print('SHA-256:', hashlib.sha256(args.output.read_bytes()).hexdigest())
    print('Run both independent certificate replays before registering a bound.')


if __name__ == '__main__':
    main()
