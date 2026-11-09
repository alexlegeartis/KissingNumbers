"""Exact-sign regression checks and semantic corrupted-witness rejection tests."""
import contextlib
import copy
import decimal
import importlib.util
import io
import json
import random
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('checker', ROOT / 'source/verify_general.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def main():
    decimal.getcontext().prec = 120
    Q = decimal.Decimal
    roots = [Q(1), Q(2).sqrt(), Q(3).sqrt(), Q(6).sqrt()]
    rng = random.Random(91823)
    for _ in range(10000):
        a = [rng.randrange(-10**20, 10**20) for _ in range(4)]
        value = sum(Q(x) * r for x, r in zip(a, roots))
        if checker.field_sign(a) != ((value > 0) - (value < 0)):
            raise AssertionError('Biquadratic sign regression')
    for a in range(-150, 151):
        for b in range(-150, 151):
            value = Q(a) + Q(b) * roots[2]
            if checker.quadratic_sign(a, b) != ((value > 0) - (value < 0)):
                raise AssertionError('Quadratic sign regression')
    print('10,000 biquadratic signs and 90,601 quadratic signs passed.', flush=True)
    for k in range(1, 101):
        p = 10**k
        for index in (1, 2, 3):
            n = int((roots[index] * Q(p)).to_integral_value(rounding=decimal.ROUND_FLOOR))
            a = [-n, 0, 0, 0]
            a[index] = p
            if checker.field_sign(a) != 1:
                raise AssertionError('Positive near cancellation')
            a[0] -= 1
            if checker.field_sign(a) != -1:
                raise AssertionError('Negative near cancellation')
    print('600 near-cancellation signs passed.', flush=True)

    original = json.loads((ROOT / 'source/certificate25.json').read_text())
    tests = []
    c = copy.deepcopy(original)
    c['points'] += 1
    tests.append(('wrong_count', c))
    c = copy.deepcopy(original)
    c['small'][0]['numerator'][0] = str(int(c['small'][0]['numerator'][0]) + 1)
    tests.append(('wrong_norm', c))
    c = copy.deepcopy(original)
    c['small'][1] = copy.deepcopy(c['small'][0])
    tests.append(('duplicated_point', c))
    c = copy.deepcopy(original)
    c['local_owners'].pop(0)
    c['points'] += 1
    tests.append(('retained_blocker', c))
    c = copy.deepcopy(original)
    c['rotation_half_tangent_numerators'].pop()
    tests.append(('wrong_rotation_count', c))
    c = copy.deepcopy(original)
    c['small'][0]['denominator'] = '0'
    tests.append(('zero_denominator', c))
    c = copy.deepcopy(original)
    c['classes'][0]['shift'][0] += 1
    tests.append(('wrong_circle_shift', c))
    with tempfile.TemporaryDirectory(prefix='kissing25_negative_') as directory:
        path = Path(directory) / 'corrupted.json'
        output = Path(directory) / 'result.json'
        for name, cert in tests:
            path.write_text(json.dumps(cert))
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    checker.verify(path, output)
            except AssertionError as error:
                print(name + ' correctly rejected: ' + str(error), flush=True)
            else:
                raise AssertionError('Corrupted certificate accepted: ' + name)
    print('All seven corrupted witnesses rejected; all self-tests passed.', flush=True)


if __name__ == '__main__':
    main()
