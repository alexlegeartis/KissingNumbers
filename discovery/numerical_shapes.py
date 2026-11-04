"""Continuous proposals that change every shape in several sign orbits.

The numerical optimizer has no certification role. It constrains every group
difference and every retained-base point, then reports the actual worst cosine.
"""
import numpy as np

from .engines.characters import _binary_rows


def orbit_shape_objective(generators, shapes, base, context):
    """Build a differentiable hinge loss and a separate full residual audit.

    Shapes use Euclidean coordinates. The base must be exact integer rows and
    closed under the stated finite binary sign action; closure is checked.
    Whole-orbit self-pairs include every nonidentity action, so free sign action
    is required and never inferred from a numerical near-duplicate test.
    """
    if not isinstance(shapes, list) or not shapes or not isinstance(shapes[0], list) or not shapes[0]:
        raise ValueError('Nonempty shape array required')
    d, count = len(shapes[0]), len(shapes)
    context.limit('max_dimension', d)
    if not isinstance(generators, list) or not generators:
        raise ValueError('Nonempty sign generators required')
    _binary_rows(generators, d, 'generators')
    size = 1 << len(generators)
    context.limit('max_codewords', size)
    terms = size*count*(count+1)//2-count
    context.limit('max_gram_entries', max(terms*d, count*len(base)))
    context.limit('max_materialized_points', count*size+len(base))
    context.limit('max_rank_work', len(base)*size*d)
    a = np.array(shapes, float)
    norms = np.linalg.norm(a,axis=1)
    if a.shape != (count, d) or not np.isfinite(a).all() or not np.isfinite(norms).all() or np.any(norms == 0):
        raise ValueError('Shapes must be finite nonzero rectangular rows')
    columns = [sum(g[j] << i for i,g in enumerate(generators)) for j in range(d)]
    signs = np.array([[1-2*((m&col).bit_count()%2) for col in columns] for m in range(size)])
    if any(np.any(np.all(a*signs[m] == a, axis=1)) for m in range(1,size)):
        raise ValueError('Continuous shape search requires free sign orbits')
    if not isinstance(base, list) or any(not isinstance(row, list) or len(row)!=d or
            any(type(x) is not int for x in row) for row in base):
        raise ValueError('Retained base requires exact integer rows')
    base_set = set(map(tuple,base))
    # Python integer signs avoid overflow in closure, even for huge base rays.
    for row in base:
        if not any(row):
            raise ValueError('Retained base has a zero vector')
        for m in range(size):
            if tuple(x*int(s) for x,s in zip(row,signs[m])) not in base_set:
                raise ValueError('Retained base is not invariant under the sign group')
    fixed = np.array(base,float).reshape((-1,d))
    if not np.isfinite(fixed).all():
        raise ValueError('Base is outside the numerical proposal range')
    if len(fixed):
        norms = np.linalg.norm(fixed,axis=1)
        if not np.isfinite(norms).all() or np.any(norms==0):
            raise ValueError('Base norm is outside the numerical proposal range')
        fixed /= norms[:,None]
    ii,jj,gg=[],[],[]
    for i in range(count):
        for j in range(i+1):
            for m in range(size):
                if i==j and m==0:
                    continue
                ii.append(i);jj.append(j);gg.append(m)
    ii,jj=np.array(ii,int),np.array(jj,int)
    sigma=signs[np.array(gg,int)]

    def residuals(flat):
        raw=np.asarray(flat,float).reshape((count,d))
        norms=np.linalg.norm(raw,axis=1)
        if not np.isfinite(norms).all() or np.any(norms==0):
            raise ValueError('Optimizer proposed zero/nonfinite shapes')
        unit=raw/norms[:,None]
        pairs=np.einsum('ij,ij,ij->i',unit[ii],unit[jj],sigma)
        base_products=unit@fixed.T
        return raw,norms,unit,pairs,base_products

    def objective(flat):
        raw,norms,unit,pairs,bp=residuals(flat)
        rp=np.maximum(pairs-.5,0.);rb=np.maximum(bp-.5,0.)
        loss=float(rp@rp+np.sum(rb*rb))
        gradient=np.zeros_like(unit)
        np.add.at(gradient,ii,2*rp[:,None]*sigma*unit[jj])
        np.add.at(gradient,jj,2*rp[:,None]*sigma*unit[ii])
        if len(fixed):
            gradient+=2*rb@fixed
        gradient=(gradient-unit*np.sum(gradient*unit,axis=1)[:,None])/norms[:,None]
        return loss,gradient.ravel()

    def audit(flat):
        raw,norms,unit,pairs,bp=residuals(flat)
        pairmax=float(np.max(pairs,initial=-1.));basemax=float(np.max(bp,initial=-1.))
        return {'max_pair_cosine':pairmax,'max_base_cosine':basemax,
                'pair_constraint_count':len(pairs),'base_constraint_count':int(bp.size),
                'shapes':unit.tolist(),'proposed_orbit_points':count*size,
                'verification_status':'numerical','new_sota':False,
                'scope':'Floating point proposal; every retained-base image checked by exact action closure. No exact spherical certificate.'}
    return objective,audit
