"""Serializable construction artifacts and bounded execution context.

A stored Verification is a receipt, never authority: verify must replay a checker.
"""
from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import math

DEFAULT_BUDGET = {
    'max_materialized_points': 20000, 'max_pairs': 10000000,
    'max_codewords': 262144, 'max_orbit_points': 10000,
    'max_trials': 200, 'solver_seconds': 30, 'max_dgls_words': 16777216,
    'max_dimension': 256, 'max_nodes': 100,
    'max_gram_entries': 2000000, 'max_rank_work': 50000000,
    'max_rational_bits': 4096,
}


def json_value(value):
    """Normalize NumPy scalars/arrays without permitting NaN or arbitrary objects."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError('Nonfinite numbers are not serializable mathematical evidence')
        return value
    if isinstance(value, dict):
        if any(not isinstance(k,str) for k in value):
            raise ValueError('Object keys must be strings')
        return {k:json_value(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [json_value(x) for x in value]
    if type(value).__module__.startswith('numpy'):
        return json_value(value.tolist() if hasattr(value,'tolist') else value.item())
    raise ValueError(f'Unsupported JSON object {type(value).__name__}')


def canonical_json(value):
    return json.dumps(json_value(value),sort_keys=True,separators=(',',':'),allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


@dataclass
class Verification:
    status: str = 'unverified'
    method: str = 'No checker has run'
    details: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.status not in ('unverified','numerical','exact'):
            raise ValueError(f'Unknown verification status {self.status!r}')
        if not isinstance(self.method,str) or not isinstance(self.details,dict):
            raise ValueError('Invalid verification fields')

    def to_dict(self):
        return json_value({'status':self.status,'method':self.method,'details':self.details})


@dataclass
class Artifact:
    kind: str
    data: dict
    facts: dict = field(default_factory=dict)
    verification: Verification = field(default_factory=Verification)

    def __post_init__(self):
        if not isinstance(self.kind,str) or not self.kind or not isinstance(self.data,dict) or not isinstance(self.facts,dict):
            raise ValueError('Artifact requires kind, data object and facts object')
        if not isinstance(self.verification,Verification):
            raise ValueError('Artifact verification must be a Verification')

    def to_dict(self):
        return json_value({'kind':self.kind,'data':self.data,'facts':self.facts,'verification':self.verification.to_dict()})

    @classmethod
    def from_dict(cls,raw):
        if not isinstance(raw,dict) or set(raw)-{'kind','data','facts','verification'} or not {'kind','data'}<=set(raw):
            raise ValueError('Invalid artifact schema')
        v=raw.get('verification',{})
        if not isinstance(v,dict) or set(v)-{'status','method','details'}:
            raise ValueError('Invalid verification schema')
        # Normalize and reject NaN even if json.load accepted it.
        raw=json_value(raw)
        return cls(raw['kind'],raw['data'],raw.get('facts',{}),Verification(**v))


@dataclass
class RunContext:
    workdir: Path
    budget: dict = field(default_factory=dict)

    def __post_init__(self):
        self.workdir=Path(self.workdir)
        merged=dict(DEFAULT_BUDGET)
        for k,v in self.budget.items():
            if k not in merged:
                raise ValueError(f'Unknown budget field {k!r}')
            if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<=0:
                raise ValueError(f'Budget {k} must be positive and finite')
            if k!='solver_seconds' and not isinstance(v,int):
                raise ValueError(f'Budget {k} must be an integer')
            merged[k]=v
        self.budget=merged

    def limit(self,key,value):
        if value>self.budget[key]:
            raise ValueError(f'{key} budget exceeded: {value} > {self.budget[key]}')


def check_params(params,allowed):
    if not isinstance(params,dict) or set(params)-set(allowed):
        raise ValueError(f'Unknown parameters: {set(params)-set(allowed) if isinstance(params,dict) else params!r}')


def integer(value,name,minimum=1):
    if isinstance(value,bool) or not isinstance(value,int) or value<minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value
