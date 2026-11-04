"""Build/check exactly24 source-backed coverage rows; no network or installations.

--write executes eligible integrated small controls, saves receipts, and writes matrix.
--check replays saved integrated proofs and validates source/evidence/recipe hashes.
Default emits the24 saved rows as JSON to stdout after checking them.
"""
import sys
sys.dont_write_bytecode=True
import argparse,hashlib,html,json,pathlib,re,time,datetime,gzip
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT))
from discovery.core import RunContext,DEFAULT_BUDGET
from discovery.runner import registry,run_program,output_artifact,validate_program
from discovery.certification import verify_artifact
BASE=HERE.parent/'baselines'
ALLOWED={'groups.ganzhinov10','groups.ganzhinov14','lattice.root_system','geometry.gram_code','groups.real_stabilizer_sector','search.support_packing'}
EXTERNAL_RECIPES=('research/coupled13/program_H2.json',
 'research/orbit_controls_20261104/control17_program.json',
 'research/orbit_controls_20261104/control18_program.json',
 'research/character_families_20261104/control19_program.json',
 'research/orbit_controls_20261104/leech_sections_program.json')
RESOURCE_PROFILES={
 'research/orbit_controls_20261104/control17_program.json':{'max_pairs':20000000},
 'research/orbit_controls_20261104/control18_program.json':{'max_pairs':40000000},
 'research/character_families_20261104/control19_program.json':{'max_pairs':80000000,'max_materialized_points':15000},
 'research/orbit_controls_20261104/leech_sections_program.json':{
  'max_materialized_points':250000,'max_orbit_points':250000,
  'max_rank_work':1000000000,'max_pairs':250000}}

def recipe_context(path):return RunContext(HERE,RESOURCE_PROFILES.get(rel(path),{}))

def read_receipt(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rb') as stream:
            raw=stream.read(200000001)
            if len(raw)>200000000:raise ValueError('Large-control receipt exceeds decompressed byte cap')
        return json.loads(raw)
    return json.loads(path.read_text())
FIXTURES={f'roots_{name}.json':d for d,name in {1:'a1',2:'a2',3:'d3',4:'d4',5:'d5',6:'e6',7:'e7',8:'e8'}.items()}
EXPECTED=[2,6,12,24,40,72,126,240,306,510,604,841,1154,1932,2564,4320,5730,8358,12460,19448,31560,49896,93150,196560]
PROTOTYPES={
 5:[('research/discovery_campaign_20261104/groups/d5_40_exact_gram.json',40,'search-discovered numerical control reconstructed as an exact Gram certificate; prototype unless a registered recipe imports it'),('research/discovery_campaign_20261104/groups/d5_40_exact_receipt.json',40,'reported independent exact Fraction Gram/rank/pair replay; source evidence only until integrated replay')],
 9:[('verifications/closed/dim09-19-record-maximality/README.md',306,'record-coordinate maximality study; not an integrated construction recipe')],
 11:[('verifications/closed/dim09-19-record-maximality/README.md',604,'record-coordinate maximality study; not an integrated construction recipe')],
 12:[('verifications/closed/dim09-19-record-maximality/README.md',841,'record-coordinate maximality study; not an integrated construction recipe')],
 13:[('verifications/closed/dim09-19-record-maximality/README.md',1154,'record-coordinate maximality study; not an integrated construction recipe')],
 15:[('verifications/closed/dim09-19-record-maximality/README.md',2564,'record-coordinate maximality study; not an integrated construction recipe')],
 16:[('research/discovery_campaign_20261104/groups/results.json',4320,'real-stabilizer-family formula/control at16D; independent representative enumeration was32D; not an integrated16D recipe'),('research/discovery_campaign_20261104/groups/README.md',4320,'family proof scope; not a separate exhaustive16D all-pairs receipt')],
 17:[('verifications/closed/dim17-layered-family/README.md',5730,'exact ceiling and matching known construction within layered family; not an integrated recipe')],
 18:[('verifications/improved/dim18-bent-hexagon/README.md',8358,'existing exact reconstruction package; separate from discovery registry'),('verifications/improved/dim18-bent-hexagon/verify18.py',8358,'exact integer plus quadratic-field pair verifier; not rerun by matrix builder')],
 19:[('verifications/closed/dim17-23-cohn-li-mechanism/README.md',11948,'historical11948 prototype/mechanism is below current12460; no current target reproduction')],
 20:[('verifications/closed/dim17-24-layer-identity/README.md',19448,'separate exact layered reconstruction, not integrated recipe')],
 21:[('verifications/closed/dim17-24-layer-identity/README.md',29768,'historical29768 reconstruction is below current31560; no current target reproduction')],
 22:[('verifications/closed/dim22-23-maximal-cross-sections/README.md',49896,'Leech section construction/obstruction package; not integrated recipe')],
 23:[('verifications/closed/dim22-23-maximal-cross-sections/README.md',93150,'Leech section construction/obstruction package; not integrated recipe')],
 24:[('common/leech.py',196560,'older Golay/Leech formula builder; the new registered exact-orbit recipe separately supplies compressed separation replay')],
}
GAPS={9:'Leech/laminated section or other306-point construction with exact coordinates/proof',11:'current604-point construction and independently replayable proof',12:'current841-point construction and independently replayable proof',13:'Zinoviev–Ericson1154-point construction and independently replayable proof',15:'2564-point laminated construction and independently replayable proof',16:'Barnes–Wall4320 target constructor plus integrated exact proof',17:'Cohn–Li5730 construction as composable program',18:'existing8358 construction adapter with exact quadratic-field proof transfer',19:'current12460 construction; historical11948 mechanism is insufficient',20:'existing19448 layered constructor adapter with exact cross-layer proofs',21:'current31560 construction; historical29768 recipe is insufficient',22:'Leech section with count49896 and compressed/exact separation proof',23:'Leech section with count93150 and compressed/exact separation proof',24:'Golay/Leech196560 constructor plus compressed separation certificate beyond generic pair budget'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
def load_source():
    meta=json.loads((BASE/'retrieval.json').read_text())['cohn'];raw=(BASE/meta['saved_path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=meta['sha256']:raise ValueError('Cohn byte hash mismatch')
    table=re.search(r'<table\b.*?</table>',raw.decode(),re.S).group()
    rows={}
    for tr in re.findall(r'<tr\b.*?</tr>',table,re.S):
        cells=re.findall(r'<t[hd][^>]*>(.*?)</t[hd]>',tr,re.S)
        text=[' '.join(html.unescape(re.sub('<[^>]+>',' ',c)).split()) for c in cells]
        if text and text[0].isdigit():
            d=int(text[0]);refs=[]
            for a in re.findall(r'<a\b[^>]*>',cells[-1]):
                title=re.search(r'title="([^"]*)"',a);href=re.search(r'href="([^"]*)"',a)
                refs.append({'citation':html.unescape(title.group(1)) if title else None,'href':meta['url']+(href.group(1) if href else '')})
            rows[d]={'value':int(text[1]),'references':refs,'reference_labels':text[-1]}
    if [rows[d]['value'] for d in range(1,25)]!=EXPECTED:raise ValueError('Source targets changed; inspect before updating expected values')
    return meta,rows

def recipe_origin(program,output_id):
    by_id={node['id']:node for node in program['nodes']};leaves=[];seen=set()
    def visit(name):
        if name in seen:return
        seen.add(name);node=by_id[name]
        if not node.get('inputs'):leaves.append(node['op'])
        for parent in node.get('inputs',[]):visit(parent)
    visit(output_id)
    imported=any(op in {'geometry.integer_code','geometry.gram_code','geometry.quadratic_code'} for op in leaves)
    generated=any(op not in {'geometry.integer_code','geometry.gram_code','geometry.quadratic_code'} for op in leaves)
    return 'mixed' if imported and generated else ('imported' if imported else 'generated')

def recipe_files():
    # Explicitly audited research programs use the same registered runner language.
    return sorted(list((ROOT/'discovery/recipes').glob('*.json'))+[ROOT/p for p in EXTERNAL_RECIPES])

def candidates():
    ops,_=registry();items=[];ctx=RunContext(HERE)
    for p in recipe_files():
        program=json.loads(p.read_text());nodes=program.get('nodes',[])
        # Audit only exact controls, not arbitrary expensive/numerical searches.
        selected=rel(p) in EXTERNAL_RECIPES or (len(nodes)==1 and nodes[0].get('op') in ALLOWED)
        selected=selected or (bool(nodes) and all(n.get('op') in {'lattice.root_system','geometry.section','geometry.gram_code'} for n in nodes))
        selected=selected or p.name in {'octacode.json','orthogonal_union.json','syndrome_cosets.json'} or p.name in FIXTURES
        if p.name in FIXTURES:
            from discovery.benchmarks.root_controls import root_recipe
            if program!=root_recipe(FIXTURES[p.name]):raise ValueError('Root fixture does not reproduce formula generator')
        if not selected:continue
        if any(n.get('op') not in ops for n in nodes):continue
        validate_program(program,ctx)
        items.append((p,program))
    return items

def write():
    meta,source=load_source();ctx=RunContext(HERE);receipts=HERE/'receipts';receipts.mkdir(exist_ok=True)
    integrated=[]
    for p,program in candidates():
        ctx=recipe_context(p)
        started=time.monotonic();receipt=run_program(program,ctx)
        receipt_path=receipts/(p.stem+'.json')
        if rel(p) in RESOURCE_PROFILES:
            receipt_path=receipt_path.with_suffix('.json.gz')
            with gzip.open(receipt_path,'wt',encoding='utf8',compresslevel=3) as stream:
                json.dump(receipt,stream,separators=(',',':'));stream.write('\n')
        else:receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
        for output_id in program['outputs']:
            artifact=output_artifact(receipt,output_id);report=verify_artifact(artifact,ctx)
            facts=report.details.get('facts',report.details)
            d=facts.get('dimension',artifact.facts.get('dimension'));count=facts.get('cardinality',artifact.facts.get('cardinality'))
            if d is None or count is None:raise ValueError(f'Missing dimension/cardinality for {p}')
            if not 1<=d<=24:raise ValueError(f'Eligible control is outside1–24; review audit scope: {p}')
            if report.status!='exact':raise ValueError(f'Eligible exact control was not exact: {p}')
            op=next(n['op'] for n in program['nodes'] if n['id']==output_id)
            origin=recipe_origin(program,output_id);provenance='generated known construction control; no autonomous discovery claim' if origin=='generated' else 'supplied exact coordinate/Gram source in reachable ancestry; geometric composition does not change its import provenance'
            if p.name in FIXTURES:provenance='literal coordinate input to executor; independently regenerated from root_controls.root_recipe formulas with exact recipe equality; no record-coordinate file or search input'
            if op=='geometry.gram_code':provenance='exact supplied Gram certificate; source numerical discovery history does not establish autonomous program synthesis'
            if op=='groups.ganzhinov10':provenance='generated known exact recipe; separate numerical subgroup-character search recovered line systems; exact recipe is reconstruction'
            if op=='groups.real_stabilizer_sector':provenance='generated implicit affine-flat/quadratic-sign sector; symbolic count and separation replay; no supplied record coordinates or novelty claim'
            if rel(p)=='research/coupled13/program_H2.json':provenance='supplied classical equator coordinates and quadratic-field cap coefficients; registered conditional switching/free-zone joining/coupled lift; known1154 construction, no novelty claim'
            if p.name=='control17_program.json':provenance='literal runtime coordinate input independently generated by build17.py from Boolean quadratic/RM code formulas; all16413585 weighted integer pairs verified; known5730 control, no held-out discovery claim'
            if p.name=='control18_program.json':provenance='native weighted quadratic coordinate input reconstructed by build18.py from the existing bent-hexagon component codes; all34923903 pairs verified; known8358 control, no held-out discovery claim'
            if p.name=='leech_sections_program.json':provenance='runtime exact orbit generation from four Golay/sextet actions; complete closure/transitivity/seed-separation proof and exact subset sections; known controls, no held-out discovery or novelty claim'
            if op=='groups.ganzhinov14':provenance='generated known E7-frame recipe; not autonomous U(3,3) subgroup-character discovery'
            if op=='search.support_packing':provenance='runtime triple-disjoint support-orbit search, with no target cardinality or record coordinates; separate blind1_24 benchmark records held-out rediscovery; geometry alone does not authenticate history'
            integrated.append({'recipe':rel(p),'recipe_sha256':sha(p),'output_id':output_id,'operators':[n['op'] for n in program['nodes']],
              'dimension':d,'ambient_dimension':facts.get('ambient_dimension',artifact.facts.get('ambient_dimension')), 'cardinality':count,'origin':origin,'provenance':provenance,
              'receipt':rel(receipt_path),'receipt_sha256':sha(receipt_path),'independent_verification':{'status':'exact_replayed','method':report.method,'scope':'current verifier rerun separately after shared-runner terminal replay; not historical-discovery or novelty proof'},
              'artifact_kind':artifact.kind,'proof_summary':{key:report.details[key] for key in ('materialized_spherical_points','independent_count_method','proof_scope','rank_scalar_work','pair_cases','pairs_checked','rank_lower_bound','arithmetic') if key in report.details},'run_seconds':round(time.monotonic()-started,6),'resource_limits':dict(ctx.budget)})
    rows=[]
    for d in range(1,25):
        target=source[d]['value'];available=[x for x in integrated if x['dimension']==d];matched=[x for x in available if x['cardinality']==target]
        external=[]
        for path,count,note in PROTOTYPES.get(d,[]):
            p=ROOT/path
            if p.exists():external.append({'path':path,'sha256':sha(p),'reported_cardinality':count,'integration':'not credited as integrated coverage','verification_status':'existing evidence inspected; not rerun by this matrix','note':note})
        coverage='exact_target_recipe' if matched else ('partial_recipe_only' if available else ('prototype_only' if external else 'none'))
        rows.append({'dimension':d,'sourced_target':target,'target_kind':'proved exact kissing number' if d in {1,2,3,4,8,24} else 'current sourced lower-bound construction; no optimality claim',
          'source':{'url':meta['url'],'retrieved_at_utc':meta['attempt_utc'],'http_date':meta.get('response_date'),'html_path':rel(BASE/meta['saved_path']),'html_sha256':meta['sha256'],'table_path':rel(BASE/'cohn_table.json'),'references':source[d]['references'],'reference_labels':source[d]['reference_labels']},
          'coverage_status':coverage,'implemented_recipe_or_none':[x['recipe'] for x in matched] or None,
          'integrated_runs':available,'implemented_numerical_search_recipes':([{'recipe':'discovery/recipes/discover10.json','sha256':sha(ROOT/'discovery/recipes/discover10.json'),'origin':'search-discovered numerical candidate','independent_verification':'not exact; bounded numerical search recipe not rerun by matrix','resource_parameters':json.loads((ROOT/'discovery/recipes/discover10.json').read_text())['nodes'][-1]['params']}] if d==10 else []),'imported_generated_search_discovered':sorted({x['origin'] for x in matched}) if matched else [],
          'independent_verification_status':'exact_target_replayed' if matched else 'no integrated exact target proof',
          'separate_repository_evidence':external,'remaining_gap':None if matched else GAPS.get(d,'current target constructor and exact proof registered in shared engine'),
          'resource_limits':{'shared_default':dict(DEFAULT_BUDGET),'target_explicit_pair_count':target*(target-1)//2,'fits_default_explicit_point_and_pair_caps':target<=DEFAULT_BUDGET['max_materialized_points'] and target*(target-1)//2<=DEFAULT_BUDGET['max_pairs'],'note':'A fitting budget does not establish an implemented recipe; compressed proof backends may avoid explicit pair limits.'},
          'new_record':False})
    output={'schema_version':1,'audit_date_context':'2026-11-04','built_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Exactly24 current Cohn targets. Integrated exact construction coverage is separate from prototypes, data import, numerical discovery and novelty. No new record.','discovery_source_sha256':{rel(p):sha(p) for p in sorted((ROOT/'discovery').rglob('*.py'))},'all_recipe_files_sha256':{rel(p):sha(p) for p in recipe_files()},'audited_eligible_recipe_paths':[rel(p) for p,_ in candidates()],'audited_eligible_outputs':[[rel(p),output] for p,program in candidates() for output in program['outputs']],'rows':rows}
    (HERE/'matrix.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
    make_readme(output)
    check();print(json.dumps({'status':'PASS','rows':24,'exact_target_dimensions':[r['dimension'] for r in rows if r['coverage_status']=='exact_target_recipe'],'prototype_or_partial_or_none_dimensions':[r['dimension'] for r in rows if r['coverage_status']!='exact_target_recipe']}))

def make_readme(matrix):
    rows=matrix['rows'];meta=rows[0]['source'];solved=[r['dimension'] for r in rows if r['target_kind'].startswith('proved')]
    header=f'''# Executable coverage matrix, dimensions1–24

The matrix has exactly24 rows. Its source is [Cohn's live table]({meta['url']}), successfully retrieved at **{meta['retrieved_at_utc']}** through the configured workspace network path. Saved HTML: `../baselines/cohn.html`; SHA-256 **{meta['html_sha256']}**. Browser502/other earlier fetch failures do not invalidate this later successful timestamped response. No new network retrieval was needed for this matrix.

Targets are existing constructions, not new results. Only dimensions{','.join(map(str,solved))} have proved exact kissing numbers in this range. The source citations, source hash/time, recipe hashes, proof receipts, provenance and resource limits are recorded per row in `matrix.json`. Separate repository prototypes receive no integrated-coverage credit. Explicitly audited research recipes that execute through the registered runner are eligible; the research allow-list includes `research/coupled13/program_H2.json` plus the17D/18D and Leech orbit/section programs under `research/orbit_controls_20261104/`, and the rational19D12460 variant under `research/character_families_20261104/`. Explicit per-recipe resource profiles govern the larger controls, whose receipts use bounded gzip serialization. Planner configuration files without runner nodes and UNVERIFIED forecasts are excluded. The existing CLI `compare` catalogue is intentionally unchanged and may remain partial (for example, a5D recipe can verify exactly while default `compare` reports no recorded5D baseline). This matrix compares replayed constructions against the captured, byte-verified source; it makes no new network retrieval.

| Dimension | Sourced target | Integrated exact target recipe | Provenance | Verification / gap |
|---|---:|---|---|---|
'''
    for r in rows:
        recipes=', '.join(pathlib.Path(p).name for p in (r['implemented_recipe_or_none'] or [])) or 'none'
        origins=', '.join(r['imported_generated_search_discovered']) or 'none'
        status='exact proof replayed' if recipes!='none' else r['coverage_status']+': '+r['remaining_gap']
        header+=f"| {r['dimension']} | {r['sourced_target']} | {recipes} | {origins} | {status} |\n"
    header+='''
`generated` means the recipe constructs a known object; it does not mean discovery of a new construction. `imported` means a certificate/coordinates are supplied as data and verified; this classification is inherited through reachable recipe ancestry. `mixed` combines supplied-coordinate/Gram and generated source ancestors. A numerical search history is recorded separately and is not retroactively turned into autonomous exact recipe discovery. The10D exact recipe is a known-result control with separate numerical line-system rediscovery. The14D exact recipe replays a supplied E7-frame construction. Where a5D Gram import is available, its exactness certifies that data, while the original randomized search provenance remains separately recorded.

Static `roots_a1` through `roots_e8` recipes are literal integer-coordinate inputs at executor level, so their ancestry is classified as imported. The audited inventory also requires exact equality with the independent `root_controls.root_recipe` formula generator. This verified generation lineage is separate from the native `lattice.root_system` routes (`roots01` through `roots08`), which construct roots during program execution. Both routes certify known constructions. The external author's original patch, source commit and author are preserved under `integration/`.

The16D octacode sign recipe produces256 points, below4320; it remains a subtarget recipe even when another16D constructor reaches the target. A registered stabilizer16 control, when present in this matrix, generates4320 points implicitly from affine-flat and quadratic-sign formulas, with symbolic count/separation proof replay rather than a4320-point array. The4D orthogonal-union demonstration imports8 supplied coordinate points, below24. The8D quotient-coset demonstration generates32 points, below240. The exact E8 section chain has target126/rank7 and72/rank6 while retaining8ambient coordinates; output IDs and ambient dimensions are recorded rather than conflated with span rank. Generic integer/Gram import and generic section capability do not imply that every target construction has been encoded. The older19D11948 and21D29768 mechanism notes are below current targets12460 and31560. The registered research `program_H2` route reaches13D1154 by conditional switching, free-zone joining and a coupled lift in Q(sqrt(3)); its supplied coordinate/cap data retain imported ancestry. Current independent replay checks all665281 pairs and exact rank13. It reproduces a known construction. Known17D5730 is reconstructed from Boolean quadratic/RM formulas and verified through the native weighted integer input. The18D8358 construction now uses the native weighted Q(sqrt(3)) backend. A rational variant of Li's known19D12460 construction uses the same native unequal-norm backend and independently checks all77619570 pairs; its printed coset choices retain imported ancestry. A four-action exact orbit recipe generates24D196560; closure, transitivity and196559 seed comparisons cover every source pair. Exact sections supply23D93150 and22D49896. The same signed-pair-normal search reaches only21D27720 and20D17400, below the captured31560 and19448 targets. Maximality of a supplied configuration neither proves that it is globally optimal nor supplies an integrated construction.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --write
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py --check
python research/discovery_campaign_20261104/coverage_1_24/build_matrix.py > /tmp/coverage_rows_1_24.json
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coverage_1_24/test_matrix.py
```

`--write` executes eligible currently registered exact controls under their declared resource profiles, separately re-verifies them, writes local receipts and rebuilds the matrix. It refuses to promote missing/unregistered/numerical recipes. `--check` verifies exactly24 ordered dimensions, unchanged source bytes, target/provenance consistency, recipe/evidence/receipt hashes, registry validity and current independent proof replay. Source fingerprints make implementation changes visible. The focused regression checks reject missing rows, a changed source target, false discovery provenance, cross-dimension attribution, omitted integrated runs and an unaudited recipe inventory. Default prints exactly24 row objects as one JSON array after checking. The larger17D/18D and Leech proofs run with explicit profiles; no installations, network requests, pushes or record-registry writes occur.

One shared construction language is used across all dimensions. Resource limits and available proof backends determine feasibility. A gap here describes missing target reproduction, not a theorem of mathematical impossibility or absence of all generic operators. Known controls and benchmark coverage remain distinct from the objective of a new exact improvement.
'''
    (HERE/'README.md').write_text(header)

def check(matrix=None):
    meta,source=load_source()
    if matrix is None:matrix=json.loads((HERE/'matrix.json').read_text())
    if matrix.get('discovery_source_sha256')!={rel(p):sha(p) for p in sorted((ROOT/'discovery').rglob('*.py'))}:raise ValueError('Discovery source fingerprints changed; rebuild matrix')
    rows=matrix['rows']
    if len(rows)!=24 or [r['dimension'] for r in rows]!=list(range(1,25)):raise ValueError('Require exactly24 rows ordered1–24')
    if matrix.get('audited_eligible_recipe_paths')!=[rel(p) for p,_ in candidates()]:raise ValueError('Eligible integrated recipe inventory changed; rebuild matrix')
    if matrix.get('all_recipe_files_sha256')!={rel(p):sha(p) for p in recipe_files()}:raise ValueError('Recipe file inventory changed; rebuild matrix')
    seen={};collected_recipes=[];ctx=RunContext(HERE)
    for r in rows:
        if r['sourced_target']!=source[r['dimension']]['value'] or r['source']['html_sha256']!=meta['sha256'] or r['source']['retrieved_at_utc']!=meta['attempt_utc'] or r['source']['url']!=meta['url'] or r['source']['references']!=source[r['dimension']]['references']:raise ValueError('Source mismatch')
        for search in r.get('implemented_numerical_search_recipes',[]):
            if sha(ROOT/search['recipe'])!=search['sha256']:raise ValueError('Numerical recipe evidence changed')
        for ext in r['separate_repository_evidence']:
            if sha(ROOT/ext['path'])!=ext['sha256']:raise ValueError('Prototype evidence changed; rebuild matrix')
        matches=[]
        for run in r['integrated_runs']:
            if run['dimension']!=r['dimension']:raise ValueError('Run dimension does not match enclosing row')
            collected_recipes.append([run['recipe'],run['output_id']])
            if sha(ROOT/run['recipe'])!=run['recipe_sha256'] or sha(ROOT/run['receipt'])!=run['receipt_sha256']:raise ValueError('Recipe/receipt hash changed')
            ctx=recipe_context(ROOT/run['recipe'])
            program=json.loads((ROOT/run['recipe']).read_text());validate_program(program,ctx)
            receipt=read_receipt(ROOT/run['receipt'])
            if receipt['program']!=program:raise ValueError('Receipt does not execute the current recipe')
            op=next(n['op'] for n in program['nodes'] if n['id']==run['output_id'])
            origin=recipe_origin(program,run['output_id'])
            if run['origin']!=origin:raise ValueError('Imported/generated provenance mismatch')
            if run['resource_limits']!=receipt['budget']:raise ValueError('Receipt resource-limit mismatch')
            if run['resource_limits']!=ctx.budget:raise ValueError('Receipt differs from declared per-recipe resource profile')
            key=(run['receipt'],run['output_id'])
            if key not in seen:
                a=output_artifact(receipt,run['output_id']);report=verify_artifact(a,ctx)
                if report.status!='exact':raise ValueError('Independent replay not exact')
                facts=report.details.get('facts',report.details);seen[key]=(facts.get('dimension',a.facts.get('dimension')),facts.get('cardinality',a.facts.get('cardinality')))
            if seen[key]!=(run['dimension'],run['cardinality']):raise ValueError('Dimension/count mismatch')
            if run['cardinality']==r['sourced_target']:matches.append(run['recipe'])
        if (matches or None)!=r['implemented_recipe_or_none']:raise ValueError('Recipe coverage mismatch')
        if bool(matches)!=(r['coverage_status']=='exact_target_recipe'):raise ValueError('Exact target coverage mislabeled')
        origins=sorted({run['origin'] for run in r['integrated_runs'] if run['recipe'] in matches})
        if r['imported_generated_search_discovered']!=origins:raise ValueError('Target provenance mismatch')
        if r['new_record']:raise ValueError('Coverage controls are not new records')
    expected=[[rel(p),output] for p,program in candidates() for output in program['outputs']]
    if matrix.get('audited_eligible_outputs')!=expected or len(collected_recipes)!=len({tuple(x) for x in collected_recipes}) or sorted(collected_recipes)!=sorted(expected):raise ValueError('Integrated runs do not cover audited eligible inventory outputs exactly once')
    return matrix

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);g=parser.add_mutually_exclusive_group();g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write:write()
    else:
        matrix=check()
        if args.check:print(json.dumps({'status':'PASS','rows':len(matrix['rows']),'source_sha256':matrix['rows'][0]['source']['html_sha256']}))
        else:print(json.dumps(matrix['rows'],ensure_ascii=False,indent=2))
