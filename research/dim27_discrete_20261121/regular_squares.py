"""Rebuild all 96 regular classes using complete native squares in the original exterior."""

from itertools import combinations, product
import json

import numpy as np
from numba import njit

import fibre_trades as ft
from bw_refill import TAILS
from fibres import native_fibres
from source import COMPLEMENT, OCTAD


@njit
def allowed_states(scores, heads_ok, choices, planes, halves, groups, line_plane, line_half):
    result = np.zeros((len(choices), 3), dtype=np.int64)
    for index, state in enumerate(choices):
        zeros, ones = np.int64(0), np.int64(0)
        valid = True
        for body, line in state:
            if not heads_ok[body, line]:
                valid = False
                break
            for j, score in enumerate(scores[body]):
                if score <= 32:
                    continue
                same_plane = planes[j] == line_plane[line]
                if score > 64 or (not same_plane and score > 48):
                    valid = False
                    break
                if same_plane:
                    bit = np.int64(1) << groups[j]
                    if halves[j] == line_half[line]:
                        ones |= bit
                    else:
                        zeros |= bit
                    if ones & zeros:
                        valid = False
                        break
            if not valid:
                break
        if valid:
            result[index] = (1, zeros, ones)
    return result


def catalogue():
    ft.initialize()
    frames, axes = native_fibres(ft.shell)
    active = {f // 32 for f in ft.owner if f >= 0}
    regular = active - {24, 134}
    assert len(regular) == 96
    keep = ~np.isin(ft.owner // 32, list(regular))
    m = ft.mids[keep]
    halves = np.array([np.count_nonzero(ft.plane[:i] == p) for i, p in enumerate(ft.plane)])
    grouped_lines = [np.flatnonzero(ft.plane == p) for p in range(3)]
    records = []
    signed = list(product((1, -1), repeat=2))
    for c, frame in enumerate(axes):
        if c in (24, 134):
            continue
        basis = frame[frame[np.arange(len(frame)), np.argmax(frame != 0, axis=1)] > 0]
        bodies, choices = [], []
        for i, j in combinations(range(8), 2):
            start = len(bodies)
            bodies.extend((s * basis[i] + t * basis[j]) // 2 for s, t in signed)
            for lines in grouped_lines:
                for k in (0, 1):
                    choices.append(np.c_[np.arange(start, start + 4), lines[[k, 1-k, 1-k, k]]])
        bodies = np.array(bodies, dtype=np.int16)
        choices = np.array(choices, dtype=np.int32)
        w = frames[c].astype(np.float32)
        mids = bodies.astype(np.float32) @ m[:, OCTAD].astype(np.float32).T
        mids += (w @ m[:, COMPLEMENT].astype(np.float32).T).max(0)
        heads = bodies @ ft.heads[:, OCTAD].T + (w @ ft.heads[:, COMPLEMENT].astype(np.float32).T).max(0)
        heads_ok = np.stack([(heads + abs(ft.heads[:, 24:] @ tail) <= 64).all(1) for tail in TAILS], axis=1)
        gates = allowed_states(mids, heads_ok, choices, ft.plane[ft.labels[keep]],
                              halves[ft.labels[keep]], ft.groups[keep], ft.plane, halves)
        for i in np.flatnonzero(gates[:, 0]):
            b, t = choices[i].T
            records.append((c, np.c_[bodies[b], TAILS[t]], gates[i, 1:]))
    return records, frames, regular


def compatibility(records, frames):
    classes = np.array([q[0] for q in records])
    labels, starts, sizes = np.unique(classes, return_index=True, return_counts=True)
    assert max(sizes) < 64
    full = np.array([(1 << int(n)) - 1 for n in sizes], dtype=np.uint64)
    rows = np.array([q[1] for q in records], dtype=np.float32)
    gates = np.array([q[2] for q in records], dtype=np.int64)
    projections = frames.reshape(-1, 16).astype(np.float32)
    cross = (projections @ projections.T).reshape(135, 32, 135, 32).max((1, 3))
    u = rows[:, :, :8].reshape(-1, 8).T
    t = rows[:, :, 8:].reshape(-1, 3).T
    values = np.arange(len(rows)) - np.repeat(starts, sizes)
    bit = np.uint64(1) << values.astype(np.uint64)
    support = []
    for k, start in enumerate(starts):
        table = []
        for i in range(start, start + sizes[k]):
            dot = rows[i, :, :8] @ u + abs(rows[i, :, 8:] @ t)
            ok = dot.reshape(4, -1, 4).max((0, 2)) + cross[classes[i], classes] <= 64
            ok &= (gates[i, 0] & gates[:, 1]) == 0
            ok &= (gates[i, 1] & gates[:, 0]) == 0
            mask = np.bitwise_or.reduceat(np.where(ok, bit, np.uint64(0)), starts)
            mask[k] = full[k]
            table.append(mask)
        support.append(np.array(table))
    return labels, full, support


def positions(mask):
    mask = int(mask)
    while mask:
        first = mask & -mask
        yield first.bit_length() - 1
        mask ^= first


def possible(initial, support, omitted, priority):
    def search(state, todo):
        while todo:
            i = todo.pop()
            compatible = np.bitwise_or.reduce(support[i][list(positions(state[i]))], axis=0)
            reduced = state & compatible
            reduced[omitted] = state[omitted]
            if not reduced.all():
                return False
            todo.extend(np.flatnonzero(reduced != state))
            state = reduced
        # A single preferred value per class supplies a cheap tentative completion.
        values = [next(positions(v)) for v in state]
        branch = None
        for i in priority:
            if i == omitted:
                continue
            for j in priority:
                if j == omitted or (int(support[i][values[i], j]) >> values[j]) & 1:
                    continue
                choices = [u for u in (i, j) if int(state[u]).bit_count() > 1]
                if not choices:
                    return False
                branch = min(choices, key=lambda u: int(state[u]).bit_count())
                break
            if branch is not None:
                break
        if branch is None:
            return True
        for k in positions(state[branch]):
            trial = state.copy()
            trial[branch] = np.uint64(1 << k)
            if search(trial, [branch]):
                return True
        return False
    return search(initial.copy(), [i for i in range(len(initial)) if i != omitted])


def main():
    records, frames, regular = catalogue()
    labels, full, support = compatibility(records, frames)
    assert len(labels) == len(regular) + 2
    priority = sorted(range(len(labels)), key=lambda i: (labels[i] in regular, i))
    for omit in range(len(labels)):
        if possible(full, support, omit, priority):
            raise ValueError(f"A packing of 97 uniform classes survives omission {labels[omit]}")
    print(json.dumps({"eligible_classes": len(labels), "candidate_square_states": len(records),
                      "excluded_omissions": len(labels), "maximum_regular_square_classes": len(regular)}, indent=2))


if __name__ == "__main__":
    main()
