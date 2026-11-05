"""Exact real quadratic-field orthogonal-projector spherical certificates."""
from fractions import Fraction
from math import isqrt
from pathlib import Path
from ..core import RunContext,Verification,digest
from ..exact_gram import bounded_fraction,check_fraction_size


def _add(x,y):return (x[0]+y[0],x[1]+y[1])
def _mul(x,y,d):return (x[0]*y[0]+d*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def _sign(x,d):
 a,b=x
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a>0 and b>0:return 1
 if a<0 and b<0:return -1
 delta=a*a-d*b*b
 return ((delta>0)-(delta<0))*((a>0)-(a<0))


def verify_artifact(artifact,context=None):
 context=context or RunContext(Path('.'));data=artifact.data
 if artifact.kind!='quadratic_projector_code' or set(data)!={'radicand','projection','selected'}:raise ValueError('invalid quadratic projector data')
 d=data['radicand'];raw=data['projection'];selected=data['selected']
 if type(d) is not int or d<2:raise ValueError('radicand must be a positive nonsquare integer')
 context.limit('max_rational_bits',d.bit_length())
 if isqrt(d)**2==d:raise ValueError('radicand must be a positive nonsquare integer')
 if not isinstance(raw,list) or not raw or any(not isinstance(row,list) or len(row)!=len(raw) for row in raw):raise ValueError('projector must be square')
 n=len(raw);context.limit('max_materialized_points',n);context.limit('max_gram_entries',n*n);context.limit('max_rank_work',n*n*n*8)
 P=[]
 for row in raw:
  R=[]
  for z in row:
   if not isinstance(z,list) or len(z)!=2:raise ValueError('projector entries need rational and radical parts')
   R.append((bounded_fraction(z[0],context,'projector coefficient'),bounded_fraction(z[1],context,'projector coefficient')))
  P.append(R)
 if any(P[i][j]!=P[j][i] for i in range(n) for j in range(i)):raise ValueError('projector not symmetric')
 for i in range(n):
  for j in range(n):
   value=(Fraction(0),Fraction(0))
   for k in range(n):
    value=_add(value,_mul(P[i][k],P[k][j],d));check_fraction_size(value[0],context);check_fraction_size(value[1],context)
   if value!=P[i][j]:raise ValueError('matrix is not an exact projector')
 trace=(sum((P[i][i][0] for i in range(n)),Fraction(0)),sum((P[i][i][1] for i in range(n)),Fraction(0)))
 if trace[1] or trace[0].denominator!=1 or trace[0]<=0:raise ValueError('invalid projector trace')
 rank=int(trace[0]);context.limit('max_dimension',rank)
 if not isinstance(selected,list) or not selected or any(not isinstance(z,list) or len(z)!=2 or type(z[0]) is not int or z[0]<0 or z[0]>=n or type(z[1]) is not int or z[1] not in (-1,1) for z in selected) or len(set(map(tuple,selected)))!=len(selected):raise ValueError('invalid selected signed indices')
 context.limit('max_materialized_points',len(selected));context.limit('max_pairs',len(selected)*(len(selected)-1)//2)
 diag=P[selected[0][0]][selected[0][0]]
 if _sign(diag,d)<=0 or any(P[i][i]!=diag for i,s in selected):raise ValueError('selected projection norms must be common and positive')
 for k,(i,s) in enumerate(selected):
  for j,t in selected[:k]:
   x=P[i][j];gap=(2*s*t*x[0]-diag[0],2*s*t*x[1]-diag[1])
   if _sign(gap,d)>0:raise ValueError('quadratic projector separation violated')
 facts={'dimension':rank,'ambient_dimension':rank,'cardinality':len(selected),'dimension_scope':'ambient projector rank; selected span may be smaller','pairs_checked':len(selected)*(len(selected)-1)//2,'max_cosine_upper_bound':[1,2]}
 return Verification('exact','symmetric idempotent over a real quadratic field, exact trace and all normalized pair inequalities',{'facts':facts,'data_sha256':digest(data),'realization':'positive semidefinite orthogonal projector in its trace dimension'})
