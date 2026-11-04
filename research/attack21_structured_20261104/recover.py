"""Rediscover linear kernels without reading the published map or labels."""
import json, time
import numpy as np
from model import C, N, P, PARITY, DATA, labels_from_holes, verify_component_geometry

def recover(seed):
    start = time.monotonic()
    rng = np.random.default_rng(seed)
    attempts = 0
    while True:
        rows = rng.integers(1, 4096, 4)
        attempts += 1
        images = (PARITY[rows[:, None] & P] * (1 << np.arange(4))[:, None]).sum(axis=0)
        if set(images) != set(range(1, 16)):
            continue
        syndromes = (PARITY[rows[:, None] & C] * (1 << np.arange(4))[:, None]).sum(axis=0)
        holes = syndromes == 0
        assert holes.sum() == 64
        labels = labels_from_holes(holes, rng)
        return dict(seed=seed, attempts=attempts, rows=list(map(int, rows)),
                    recovery_seconds=time.monotonic()-start, holes=64, selected=960), holes, labels

if __name__ == '__main__':
    DATA.mkdir(exist_ok=True)
    results, all_holes, all_labels = [], [], []
    for seed in range(32):
        report, holes, labels = recover(seed)
        report['geometry'] = verify_component_geometry(labels)
        results.append(report); all_holes.append(holes); all_labels.append(labels)
    np.savez_compressed(DATA / 'recovery.npz', holes=all_holes, labels=all_labels)
    output = dict(tests=32, successes=32, provided_published_map=False,
                  provided_published_labels=False, assumed_structure='codimension-four linear kernels',
                  maximum_seconds=max(r['recovery_seconds'] for r in results),
                  total_attempts=sum(r['attempts'] for r in results), runs=results)
    (DATA / 'recovery.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='runs'}), flush=True)
