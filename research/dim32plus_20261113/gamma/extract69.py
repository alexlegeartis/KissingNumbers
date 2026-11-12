"""Extract an exact dual in independent orthonormal D4 coordinates."""
from pathlib import Path
import json, math
from fractions import Fraction as F
import sympy as s
import numpy as np
import dk_affine as d

P=Path(__file__).parent
J=json.loads((P/'independent69_model.json').read_text())
M=d.build(4)
value,result,obj=d.solve_obj(M,J['offsets'])
A=s.Matrix(M['rows']);b=s.Matrix(M['rhs']);c=s.Matrix(list(map(int,obj)))
active=[i for i,v in enumerate(result.x*M['colscale']) if v>1e-12]
AT=A[:,active].T
sol=s.linsolve((AT,c[active,0]))
y=list(next(iter(sol)))
symbols=set().union(*(x.free_symbols for x in y))
y=s.Matrix([x.subs({z:0 for z in symbols}) for x in y])
slack=c-A.T*y
assert min(slack)>=0,list(slack)
bound=(b.T*y)[0]
assert bound==743522304,bound
den=math.lcm(*(int(x.q) for x in y))
coef=[int(x*den) for x in y]
C=dict(schema='gamma72-d4-head-compression-v1',dimension=69,population=d.N,
       prior_bound=689452922,bound=int(bound),profiles=J['offsets'],heights=J['height'],
       moment_powers=M['tags'],dual_denominator=den,dual_multipliers=coef,
       rows=M['rows'],rhs=M['rhs'],objective=list(map(int,obj)),
       groups=[sorted(gg) for gg in M['groups']],
       proof_scope='Exact existence lower bound conditional only on the cited Gamma72 lattice minimum, population and 11-design facts.')
(P/'independent69_certificate.json').write_text(json.dumps(C,indent=2))
print(json.dumps(dict(bound=int(bound),denominator=den,multipliers=coef,
     slack=[str(x) for x in slack],active_columns=active),indent=2))
