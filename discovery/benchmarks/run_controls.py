"""End-to-end exact controls; receipts go to an explicitly chosen new directory."""
import argparse
import json
from pathlib import Path
import resource
import time
from discovery.core import RunContext
from discovery.runner import run_program,output_artifact

EXPECTED={'ganzhinov10':(10,510),'ganzhinov14':(14,1932),
          'known32':(32,352800),'dgls96':(96,2**37)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',required=True)
    parser.add_argument('--controls',nargs='+',choices=list(EXPECTED),default=list(EXPECTED))
    args=parser.parse_args();out=Path(args.outdir)
    if out.exists():parser.error('Use a new output directory; benchmark receipts are not overwritten')
    out.mkdir(parents=True)
    context=RunContext(out);rows=[]
    for name in args.controls:
        recipe=json.loads((Path(__file__).parents[1]/'recipes'/f'{name}.json').read_text())
        start=time.monotonic()
        receipt=run_program(recipe,context)
        artifact=output_artifact(receipt)
        got=artifact.facts['dimension'],artifact.facts['cardinality']
        if artifact.verification.status!='exact' or got!=EXPECTED[name]:
            raise RuntimeError(f'{name} control failed: {got}')
        row={'control':name,'dimension':got[0],'cardinality':got[1],'status':'exact',
             'seconds':round(time.monotonic()-start,3),
             'cumulative_process_peak_rss_kib_linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             'kind':artifact.kind,'receipt':name+'.json','new_sota':False}
        (out/(name+'.json')).write_text(json.dumps(receipt,indent=2)+'\n')
        rows.append(row);print(json.dumps(row),flush=True)
    summary={'scope':'known controls only; runtime observed on this executor, not a performance guarantee','controls':rows}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
