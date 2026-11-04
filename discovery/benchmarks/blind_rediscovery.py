"""A separate target-blind construction run, followed by held-out evaluation.

The search process cannot read the comparison file: it is provided only to the
separate evaluate subcommand. No record recipes or coordinates are loaded.
"""
import argparse
from fractions import Fraction
import hashlib
import html
import json
from pathlib import Path
import re
from time import monotonic

from discovery.core import RunContext, digest
from discovery.runner import run_program, output_artifact, source_fingerprints
from discovery.certification import verify_artifact


def partitions(n, minimum=2, maximum=4):
    if n == 0:
        yield []
    for size in range(minimum, min(maximum,n)+1):
        for tail in partitions(n-size, size, maximum):
            yield [size]+tail


def cycle_catalogue(n):
    choices = [[n]]
    if n > 1:
        choices.append([n-1,1])
    choices.extend(list(partitions(n))[:6])
    return [list(value) for value in dict.fromkeys(tuple(x) for x in choices)]


def search(out, seconds):
    out.mkdir(parents=True, exist_ok=False)
    ctx=RunContext(out, {'max_pairs':40000000, 'solver_seconds':max(seconds,1)})
    rows=[]
    # This catalogue is independent of any external target or incumbent count.
    for dimension in range(1,25):
        alternatives=[]
        for family in ('A','D','E'):
            if family=='D' and dimension<2 or family=='E' and dimension not in (6,7,8):
                continue  # mathematical existence conditions, not quality gates
            alternatives.append(('root_formula', 'lattice.root_system', {'family':family,'rank':dimension}))
        if dimension>=4:
            alternatives.extend(('support_search','search.support_packing',
                                 {'dimension':dimension,'cycles':cycles,'seconds':seconds})
                                for cycles in cycle_catalogue(dimension))
        if dimension>=2 and dimension&(dimension-1)==0:
            alternatives.extend(('stabilizer_formula','groups.real_stabilizer_sector',
                                 {'m':dimension.bit_length()-1,'parity':p}) for p in (0,1))
        results=[]
        for index,(origin,op,params) in enumerate(alternatives):
            program={'schema_version':1,'nodes':[{'id':'sphere','op':op,'params':params}],'outputs':['sphere']}
            start=monotonic()
            receipt=run_program(program,ctx)
            artifact=output_artifact(receipt)
            report=verify_artifact(artifact,ctx)
            facts=report.details.get('facts',report.details)
            if report.status!='exact' or facts['dimension']!=dimension:
                raise ValueError('Candidate failed exact dimension/geometry replay')
            path=out/f'd{dimension:02d}-{index:02d}.json'
            path.write_text(json.dumps(receipt,separators=(',',':'))+'\n')
            results.append({'origin':origin,'op':op,'params':params,'cardinality':facts['cardinality'],
                            'receipt':path.name,'receipt_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                            'seconds':monotonic()-start})
        rows.append({'dimension':dimension,'best':max(results,key=lambda r:r['cardinality']),
                     'best_search':max((r for r in results if r['origin']=='support_search'),
                                       key=lambda r:r['cardinality'],default=None),'attempts':results})
        summary={'schema_version':1,'phase':'target_blind_search','source_sha256':source_fingerprints(),
                 'comparison_inputs':[], 'supplied_record_coordinates':False,
                 'grammar':'A/D/E roots; searched triple-disjoint weight-four support orbits; real stabilizer parity sectors',
                 'formula_controls_are_discovery':False,'budget':ctx.budget,'rows':rows}
        (out/'search.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps({'dimension':dimension,'best':rows[-1]['best']['cardinality'],
                          'best_support_search':rows[-1]['best_search']['cardinality'] if rows[-1]['best_search'] else None}),flush=True)


def read_table(path):
    raw=path.read_bytes()
    table=re.search(r'<table\b.*?</table>',raw.decode(),re.S)
    if table is None:
        raise ValueError('Missing comparison table')
    result={}
    for row in re.findall(r'<tr\b.*?</tr>',table.group(),re.S):
        cells=re.findall(r'<t[hd]([^>]*)>(.*?)</t[hd]>',row,re.S)
        values=[' '.join(html.unescape(re.sub('<[^>]+>',' ',c)).split()) for _,c in cells]
        if not values or not values[0].isdigit():
            continue
        dimension=int(values[0]); lower=int(values[1])
        # Exact dimensions use one colspan=2 cell. Treating the following
        # ratio cell as an upper bound would corrupt neighboring growth rates.
        exact=bool(re.search(r'colspan\s*=\s*["\x27]?2',cells[1][0]))
        upper=lower if exact else int(values[2])
        if lower>upper or dimension in result:
            raise ValueError('Invalid or duplicate comparison row')
        result[dimension]={'lower':lower,'upper':upper,'proved_exact':exact}
    if sorted(result)!=list(range(1,49)):
        raise ValueError('Require all 48 current neighboring comparison rows')
    return result


def fraction(value):
    return {'numerator':value.numerator,'denominator':value.denominator,'decimal':float(value)}


def growth_diagnostic(table):
    """Gap-expansion ratios, using exact rational neighboring inputs.

    Score > 1 means upper/lower widens. It does not predict an improvement.
    The centered score distinguishes a one-step plateau from a two-step one.
    """
    rows=[]
    for n in range(2,49):
        before,current=table[n-1],table[n]
        lower=Fraction(current['lower'],before['lower'])
        upper=Fraction(current['upper'],before['upper'])
        row={'dimension':n,'neighbors':{str(n-1):before,str(n):current},
             'lower_growth':fraction(lower),'upper_growth':fraction(upper),
             'gap_expansion':fraction(upper/lower),'new_record_evidence':False}
        if n<48:
            after=table[n+1]
            row['neighbors'][str(n+1)]=after
            row['centered_two_step_gap_expansion']=fraction(
                Fraction(after['upper']*before['lower'],before['upper']*after['lower']))
        rows.append(row)
    return sorted(rows,key=lambda r:r['gap_expansion']['decimal'],reverse=True)


def evaluate(directory, source, metadata):
    summary_path=directory/'search.json'
    summary=json.loads(summary_path.read_text())
    if summary['phase']!='target_blind_search' or summary['comparison_inputs'] or summary['supplied_record_coordinates']:
        raise ValueError('Search was not declared target blind')
    table=read_table(source)
    ctx=RunContext(directory,summary['budget'])
    rows=[]
    for row in summary['rows']:
        for item in row['attempts']:
            path=directory/item['receipt']
            if hashlib.sha256(path.read_bytes()).hexdigest()!=item['receipt_sha256']:
                raise ValueError('Changed benchmark receipt')
            artifact=output_artifact(json.loads(path.read_text()))
            report=verify_artifact(artifact,ctx)
            facts=report.details.get('facts',report.details)
            if report.status!='exact' or facts['dimension']!=row['dimension'] or facts['cardinality']!=item['cardinality']:
                raise ValueError('Benchmark exact replay differs')
        n=row['dimension'];target=table[n]['lower'];best=row['best']['cardinality']
        searched=row['best_search']['cardinality'] if row['best_search'] else None
        rows.append({'dimension':n,'target':target,'best_exact':best,
                     'target_matched':best==target,'ratio':fraction(Fraction(best,target)),
                     'best_support_search':searched,'support_search_target_matched':searched==target,
                     'best_origin':row['best']['origin'],'new_record':False})
    if [r['dimension'] for r in rows]!=list(range(1,25)):
        raise ValueError('Benchmark must contain exactly ordered dimensions 1 through 24')
    meta=json.loads(metadata.read_text()) if metadata else {}
    if 'cohn' in meta:
        meta=meta['cohn']
    actual_hash=hashlib.sha256(source.read_bytes()).hexdigest()
    if meta.get('sha256') not in (None,actual_hash):
        raise ValueError('Source hash does not match metadata')
    report={'schema_version':1,'phase':'held_out_evaluation','search_sha256':hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            'comparison_source':{'path':str(source),'sha256':actual_hash,'metadata':meta},
            'target_matches':[r['dimension'] for r in rows if r['target_matched']],
            'actual_support_search_target_matches':[r['dimension'] for r in rows if r['support_search_target_matched']],
            'rows':rows,'growth_diagnostic':growth_diagnostic(table),
            'scope':'Bounded grammar benchmark; formula controls, search rediscovery and new records are distinct; no new record'}
    (directory/'evaluation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'target_matches':report['target_matches'],
                      'support_search_matches':report['actual_support_search_target_matches'],
                      'growth_priority_first10':[r['dimension'] for r in report['growth_diagnostic'][:10]]}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    run=sub.add_parser('search');run.add_argument('--out',type=Path,required=True);run.add_argument('--seconds',type=float,default=.5)
    check=sub.add_parser('evaluate');check.add_argument('directory',type=Path);check.add_argument('--source',type=Path,required=True)
    check.add_argument('--metadata',type=Path)
    args=parser.parse_args()
    if args.command=='search':search(args.out,args.seconds)
    else:evaluate(args.directory,args.source,args.metadata)


if __name__=='__main__':
    main()
