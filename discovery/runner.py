"""Versioned, allow-listed construction DAG execution. Recipes are data, not code."""
from dataclasses import dataclass
from pathlib import Path
import importlib
import json
import platform
import re
import time
from .core import Artifact,RunContext,digest,json_value
from .certification import verify_artifact

MODULES=('groups','rings','geometry','cosets','layers','lattices','gluing','stabilizers','selection','signed_layers','coset_search','binary_translate_search','block_isometry_search','signed_exchange','ring_orbit_search','group_seed_search','coset_group_search','mixed_group_search','partial_group_assemblies','quadratic_group_assemblies','cavity_repacking','construction_search','exact_orbits','characters','body_heads','partial_orbits','contact_reflections','finite_field_search','spectral_graph_search')


def input_kinds_match(inputs,expected):
    """Legacy exact kind strings or explicit finite alternatives per position."""
    if len(inputs)!=len(expected):return False
    for artifact,allowed in zip(inputs,expected):
        if isinstance(allowed,str):
            if artifact.kind!=allowed:return False
        elif isinstance(allowed,(tuple,list)) and allowed and all(isinstance(kind,str) for kind in allowed):
            if artifact.kind not in allowed:return False
        else:
            raise ValueError('Invalid registered input-kind contract')
    return True


def registry():
    ops={};kinds={}
    for name in MODULES:
        module=importlib.import_module(f'discovery.engines.{name}')
        for op,handler in module.OPERATORS.items():
            if op in ops: raise ValueError('Duplicate registered operator')
            ops[op]=handler;kinds[op]=module.INPUT_KINDS[op]
    return ops,kinds


def validate_program(program,context):
    if not isinstance(program,dict) or set(program)-{'schema_version','name','description','nodes','outputs','provenance'}:
        raise ValueError('Unknown program fields')
    if type(program.get('schema_version')) is not int or program.get('schema_version')!=1: raise ValueError('Only schema_version 1 is supported')
    nodes=program.get('nodes');outputs=program.get('outputs')
    if not isinstance(nodes,list) or not nodes: raise ValueError('Program needs nodes')
    context.limit('max_nodes',len(nodes))
    if not isinstance(outputs,list) or not outputs or any(not isinstance(x,str) for x in outputs) or len(set(outputs))!=len(outputs):
        raise ValueError('Program needs distinct output ids')
    ops,kinds=registry();by_id={}
    for node in nodes:
        if not isinstance(node,dict) or set(node)-{'id','op','inputs','params'} or not {'id','op'}<=set(node):
            raise ValueError('Invalid node fields')
        name=node['id']
        if not isinstance(name,str) or re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*',name) is None:
            raise ValueError('Node id must be a simple identifier')
        if name in by_id: raise ValueError(f'Duplicate node id {name}')
        if node['op'] not in ops: raise ValueError(f'Unsupported operator {node["op"]!r}')
        ins=node.get('inputs',[])
        if not isinstance(ins,list) or any(not isinstance(x,str) for x in ins):
            raise ValueError('Inputs must be node ids')
        if len(ins)!=len(kinds[node['op']]): raise ValueError(f'Wrong input count for {node["op"]}')
        if not isinstance(node.get('params',{}),dict): raise ValueError('Node params must be an object')
        by_id[name]=node
    if any(name not in by_id for name in outputs): raise ValueError('Missing output node')
    if any(dep not in by_id for node in nodes for dep in node.get('inputs',[])):
        raise ValueError('Missing dependency')
    ordered=[];pending=list(nodes);seen=set()
    while pending:
        ready=[node for node in pending if all(dep in seen for dep in node.get('inputs',[]))]
        if not ready: raise ValueError('Construction graph has a cycle')
        for node in ready:
            ordered.append(node);seen.add(node['id']);pending.remove(node)
    json_value(program)
    return ordered,ops,kinds


def source_fingerprints():
    root=Path(__file__).parent
    paths=list(root.glob('*.py'))+list((root/'engines').glob('*.py'))
    import hashlib
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def run_program(program,context,on_progress=None):
    nodes,ops,kinds=validate_program(program,context)
    artifacts={};events=[]
    for node in nodes:
        inputs=[artifacts[x] for x in node.get('inputs',[])]
        expected=kinds[node['op']]
        if not input_kinds_match(inputs,expected):
            raise ValueError(f'{node["op"]} needs {expected}, got {tuple(a.kind for a in inputs)}')
        if on_progress:on_progress({'node':node['id'],'op':node['op'],'state':'started'})
        start=time.monotonic()
        artifact=ops[node['op']](inputs,node.get('params',{}),context)
        if not isinstance(artifact,Artifact): raise ValueError('Operator returned an invalid artifact')
        # Force JSON normalization now rather than after an expensive whole run.
        artifact=Artifact.from_dict(artifact.to_dict())
        artifacts[node['id']]=artifact
        event={'node':node['id'],'op':node['op'],'state':'completed',
               'seconds':round(time.monotonic()-start,6),'kind':artifact.kind,
               'verification_status':artifact.verification.status,'artifact_sha256':digest(artifact.to_dict())}
        events.append(event)
        if on_progress:on_progress(event)
    # Rerun the independent verifier for exact terminal results. A serialized
    # status or inherited proof from an incompatible operation is insufficient.
    for name in program['outputs']:
        artifact=artifacts[name]
        if artifact.verification.status=='exact':
            report=verify_artifact(artifact,context)
            if report.status!='exact': raise ValueError('Exact terminal receipt could not be independently replayed')
            artifact.verification=report
    for event in events:
        event['artifact_sha256']=digest(artifacts[event['node']].to_dict())
    return {'schema_version':1,'program':program,'program_sha256':digest(program),
            'budget':context.budget,'python':platform.python_version(),'source_sha256':source_fingerprints(),
            'events':events,'artifacts':{k:v.to_dict() for k,v in artifacts.items()},
            'artifacts_sha256':{k:digest(v.to_dict()) for k,v in artifacts.items()},'outputs':program['outputs'],
            'scope':'construction experiment; no automatic record registration or worldwide-priority claim'}


def output_artifact(receipt,name=None):
    if not isinstance(receipt,dict): raise ValueError('Receipt must be an object')
    if 'kind' in receipt: return Artifact.from_dict(receipt)
    if receipt.get('schema_version')!=1 or not isinstance(receipt.get('artifacts'),dict):
        raise ValueError('Invalid construction receipt')
    if receipt.get('program_sha256')!=digest(receipt.get('program')):
        raise ValueError('Receipt program hash mismatch')
    hashes=receipt.get('artifacts_sha256')
    if not isinstance(hashes,dict) or set(hashes)!=set(receipt['artifacts']):
        raise ValueError('Receipt artifact hash manifest missing or incomplete')
    if any(hashes[k]!=digest(v) for k,v in receipt['artifacts'].items()):
        raise ValueError('Receipt artifact hash mismatch')
    for event in receipt.get('events',[]):
        if event.get('artifact_sha256')!=hashes.get(event.get('node')):
            raise ValueError('Receipt event hash mismatch')
    outputs=receipt.get('outputs',[])
    if outputs!=receipt.get('program',{}).get('outputs'):
        raise ValueError('Receipt outputs differ from program')
    if name is None:
        if len(outputs)!=1: raise ValueError('Select an output id with --output')
        name=outputs[0]
    if name not in outputs: raise ValueError('Selected id is not a terminal output')
    return Artifact.from_dict(receipt['artifacts'][name])

