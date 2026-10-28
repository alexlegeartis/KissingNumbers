"""Check that both representations reject duplicated points and a bad pure angle."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    scripts = ['partial43/verify.py', 'continuation43/verify_coordinates.py']
    report = {'status':'PASS', 'rejected_cases':{}}
    with tempfile.TemporaryDirectory(prefix='kissing43_negative_') as directory:
        root = Path(directory)
        for name in ('partial43', 'continuation43'):
            shutil.copytree(HERE.parent/name, root/name,
                            ignore=shutil.ignore_patterns('__pycache__'))
        with np.load(root/'partial43/explicit_points_43.npz') as source:
            valid = dict(source)
        with np.load(root/'continuation43/delta.npz') as source:
            delta = dict(source)
        for case in ('duplicate', 'pure_angle'):
            data = {k:v.copy() for k,v in valid.items()}
            change = {k:v.copy() for k,v in delta.items()}
            if case == 'duplicate':
                data['body_indices'][-1] = data['body_indices'][0]
                data['head_numerators'][-1] = data['head_numerators'][0]
                change['added_body_indices'][-1] = data['body_indices'][0]
                change['added_head_numerators'][-1] = data['head_numerators'][0]
            else:
                with np.load(root/'partial43/configuration_43.npz') as model:
                    direction = model['pure_directions'][0].astype(float)
                direction /= np.linalg.norm(direction)
                bad_head = np.rint(100*int(data['head_denominator'])*direction).astype(np.int64)
                data['head_numerators'][-1] = bad_head
                change['added_head_numerators'][-1] = bad_head
            np.savez_compressed(root/'partial43/explicit_points_43.npz', **data)
            np.savez_compressed(root/'continuation43/delta.npz', **change)
            rejected = []
            for script in scripts:
                run = subprocess.run([sys.executable, str(root/script)], cwd=root,
                                     capture_output=True, text=True)
                if run.returncode == 0 or 'AssertionError' not in run.stderr:
                    raise RuntimeError(f'{script} did not reject {case} as a bad certificate')
                rejected.append(script)
            report['rejected_cases'][case] = rejected
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
