"""Run with python -m discovery. No installation, network or publication side effects."""
import argparse
import json
from pathlib import Path
import sys
from .core import RunContext,canonical_json
from .runner import registry,run_program,output_artifact
from .certification import verify_artifact


def read(path):
    with open(path) as stream:return json.load(stream)


def context(args):
    limits={}
    for pair in getattr(args,'budget',[]) or []:
        if '=' not in pair: raise ValueError('Budget must be KEY=VALUE')
        key,value=pair.split('=',1)
        if key in limits:raise ValueError(f'Duplicate budget {key}')
        limits[key]=float(value) if key=='solver_seconds' else int(value)
    return RunContext(Path(getattr(args,'workdir','.') or '.'),limits)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('operators',help='list implemented operators, not roadmap items')
    run=sub.add_parser('run',help='run a versioned construction program')
    run.add_argument('recipe');run.add_argument('--out',required=True,help='new JSON receipt; existing files are never overwritten')
    run.add_argument('--budget',action='append',default=[]);run.add_argument('--workdir',default='.')
    search=sub.add_parser('search',help='enumerate an explicit finite construction grammar under resource budgets')
    search.add_argument('grammar');search.add_argument('--out-dir',required=True,help='new directory for search summary and normal run receipts')
    search.add_argument('--max-programs',type=int,default=100);search.add_argument('--start-index',type=int,default=0)
    search.add_argument('--time-limit-seconds',type=float,default=60);search.add_argument('--per-program-seconds',type=float,default=30)
    search.add_argument('--max-grammar-choices',type=int,default=100000);search.add_argument('--max-receipt-bytes',type=int,default=16777216)
    search.add_argument('--budget',action='append',default=[]);search.add_argument('--workdir',default='.')
    verify=sub.add_parser('verify',help='independently reverify an artifact or terminal receipt output')
    verify.add_argument('receipt');verify.add_argument('--output');verify.add_argument('--budget',action='append',default=[])
    compare=sub.add_parser('compare',help='compare independently verified count against dated benchmark snapshot')
    compare.add_argument('receipt');compare.add_argument('--output');compare.add_argument('--budget',action='append',default=[])
    compare.add_argument('--catalog',default=str(Path(__file__).parent/'benchmarks/catalog.json'))
    args=parser.parse_args(argv)
    try:
        if args.command=='operators':
            ops,kinds=registry();print(json.dumps([{'op':op,'inputs':list(kinds[op]),'implemented':True} for op in sorted(ops)],indent=2));return 0
        ctx=context(args)
        if args.command=='search':
            from .scheduler import search_grammar
            def progress(event):print(json.dumps(event),file=sys.stderr,flush=True)
            result=search_grammar(read(args.grammar),ctx,args.out_dir,
                max_programs=args.max_programs,start_index=args.start_index,
                time_limit_seconds=args.time_limit_seconds,per_program_seconds=args.per_program_seconds,
                max_grammar_choices=args.max_grammar_choices,max_receipt_bytes=args.max_receipt_bytes,
                on_progress=progress)
            print(json.dumps({'summary':str(Path(args.out_dir)/'search.json'),'status':result['status'],
                'attempted_programs':result['attempted_programs'],'declared_programs':result['declared_programs'],
                'next_index':result['next_index'],'best_exact_by_dimension':result['best_exact_by_dimension']}))
            return 0 if result['status']=='complete' else 2
        if args.command=='run':
            out=Path(args.out)
            if out.exists():raise ValueError(f'Refusing to overwrite existing receipt {out}')
            program=read(args.recipe)
            def progress(event): print(json.dumps(event),file=sys.stderr,flush=True)
            receipt=run_program(program,ctx,progress)
            out.parent.mkdir(parents=True,exist_ok=True)
            # Exclusive creation also handles a race with another process.
            with out.open('x') as stream:stream.write(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
            print(json.dumps({'receipt':str(out),'outputs':receipt['outputs'],
                              'statuses':{x:receipt['artifacts'][x]['verification']['status'] for x in receipt['outputs']}}));return 0
        artifact=output_artifact(read(args.receipt),args.output)
        report=verify_artifact(artifact,ctx)
        if args.command=='verify':print(json.dumps(report.to_dict(),indent=2));return 0 if report.status=='exact' else 2
        if report.status!='exact':raise ValueError('Only independently exact-certified configurations may be compared')
        if artifact.kind not in ('spherical_code','implicit_spherical_code'):
            raise ValueError('Compare requires a spherical configuration; use geometry.sign_embed for a binary code')
        facts=report.details.get('facts',report.details);d=facts.get('dimension');count=facts.get('cardinality')
        catalog=read(args.catalog);rows=catalog.get('comparisons',[])
        rows=[r for r in rows if r['dimension']==d]
        if not rows: raise ValueError(f'No recorded baseline for dimension {d}')
        best=max(rows,key=lambda r:r['cardinality'])
        print(json.dumps({'dimension':d,'cardinality':count,'recorded_baseline':best,
                          'difference':count-best['cardinality'],'exceeds_recorded_baseline':count>best['cardinality'],
                          'scope':'dated comparison only; fresh literature review and reviewable certificate required before any new-SOTA claim'},indent=2))
        return 0
    except (ValueError,KeyError,TypeError,OSError,ArithmeticError) as exc:
        print(f'discovery: {exc}',file=sys.stderr);return 1

if __name__=='__main__':raise SystemExit(main())
