"""Extract rational duals from exploratory LP active columns; verify exactly."""
from pathlib import Path
from fractions import Fraction
import json
import sys
import sympy as sp

for arg in sys.argv[1:]:
    path = Path(arg)
    d = json.loads(path.read_text())
    A = sp.Matrix(d['rows'])
    b = sp.Matrix(d['rhs'])
    c = sp.Matrix([sp.Rational(s) for s in d['objective']])
    active = [i for i, x in enumerate(d['primal']) if x > 1e-12]
    sol = next(iter(sp.linsolve((A[:, active].T, c[active, :]))))
    free = set().union(*(x.free_symbols for x in sol))
    if free:
        # Setting free multipliers to zero is a candidate, never a proof by itself.
        sol = [x.subs({f: 0 for f in free}) for x in sol]
    y = sp.Matrix(sol)
    slack = c - A.T*y
    assert all(x >= 0 for x in slack), list(slack)
    value = (b.T*y)[0]
    ceiling = int(sp.ceiling(value))
    d.update(dual_exact=[str(x) for x in y], slack_exact=[str(x) for x in slack],
             exact_value=str(value), exact_bound=ceiling, active_columns=active)
    if d.get('dimension') == 69:
        d['with_two_poles'] = ceiling + 2
    out = path.with_name(path.stem+'_exact.json')
    out.write_text(json.dumps(d, indent=2)+'\n')
    print(out.name, 'rank', A.rank(), 'active', len(active), 'value', value,
          'ceiling', ceiling, 'poles', d.get('with_two_poles'), flush=True)
