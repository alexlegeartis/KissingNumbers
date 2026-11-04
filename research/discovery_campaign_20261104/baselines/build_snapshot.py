"""Read-only comparison builder. Writes only alongside this script; no live fetches."""
import pathlib,json,re,hashlib,datetime,ast
HERE=pathlib.Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def read(p):return p.read_text()
def num(s):return int(re.sub(r'[^0-9]','',s))
rows=[]
for line in read(ROOT/'RESULTS.md').splitlines():
    cells=[x.strip() for x in line.split('|')[1:-1]]
    if len(cells)>=6 and cells[0].isdigit() and int(cells[0])>=25:
        rows.append({'dimension':int(cells[0]),'historical_comparison':num(cells[1]),'repository_count':num(cells[2]),'mechanism':cells[4].strip('`'),'construction':cells[5]})
assert len(rows)==69
cohn={int(r[0]):int(r[1]) for r in json.loads(read(HERE/'cohn_table.json'))[1:]}
# Compare fresh 5–16D records with the repository's explicitly pinned COHN table.
module=ast.parse(read(ROOT/'common/published.py'))
pinned=next(ast.literal_eval(n.value) for n in module.body
            if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='COHN' for t in n.targets))
table_rows={int(r[0]):r for r in json.loads(read(HERE/'cohn_table.json'))[1:]}
retrieval_records=json.loads(read(HERE/'retrieval.json'))
for record in retrieval_records.values():
    if record.get('status')==200:
        assert hashlib.sha256((HERE/record['saved_path']).read_bytes()).hexdigest()==record['sha256']
retrieval=retrieval_records['cohn']
small=[]
for d in range(5,17):
    cells=table_rows[d]
    solved=(d==8)
    upper=cohn[d] if solved else int(cells[2])
    assert cohn[d]==pinned[d]
    small.append({'dimension':d,'research_threshold_to_beat':cohn[d],
        'fresh_cohn_lower':cohn[d],'fresh_cohn_upper':upper,
        'kissing_number_known_exactly':solved,
        'repository_pinned_lower':pinned[d],
        'repository_registered_improvement':False,
        'source_url':retrieval['url'],'source_retrieved_at_utc':retrieval['attempt_utc'],
        'source_html_sha256':retrieval['sha256'],'source_saved_path':retrieval['saved_path'],'cohn_reference_labels':cells[-1],
        'campaign_status':'equally active shared-operator target; 8D is a solved-value control' if solved else 'equally active shared-operator research target'})
# Keep named direct competitors separate from the current Cohn entries credited to us.
competitors={25:(197583,'Takhanov-Yun'),26:(201134,'PackingStar'),27:(209760,'PackingStar'),28:(221864,'PackingStar'),29:(231676,'PackingStar'),30:(247967,'PackingStar'),31:(251169,'PackingStar'),32:(352320,'Lin'),33:(370240,'Lin'),34:(390980,'Lin'),35:(415308,'Lin'),36:(484760,'Yang et al.'),37:(500008,'Olivier'),38:(605322,'Lin'),39:(768916,'Lin'),40:(1092000,'Dorofeev-Sun-Wang'),41:(1324472,'Dorofeev-Sun-Wang'),42:(1792386,'Dorofeev-Sun-Wang'),43:(2553792,'Yang et al.'),44:(4153600,'Conway-Sloane / Dorofeev-Sun-Wang'),45:(7380090,'Yang et al.'),46:(12309600,'Conway-Sloane / Dorofeev-Sun-Wang')}
for r in rows:
    d=r['dimension'];r['research_threshold_to_beat']=r['repository_count']
    r['fresh_cohn_lower']=cohn.get(d)
    if d in cohn:assert cohn[d]==r['repository_count']
    if d in competitors:
        v,who=competitors[d];r.update(named_competitor_count=v,named_competitor=who,margin_over_named_competitor=r['repository_count']-v)
    if 49<=d<=63:r['fresh_derived_direct_sum_floor']=52416000+cohn[d-48]
    if 73<=d<=91:r['fresh_derived_direct_sum_floor']=6218175600+cohn[d-72]
    if 92<=d<=96:r['published_DGLS_shortening_floor']=2**(d-59)
    if 'fresh_derived_direct_sum_floor' in r:assert r['repository_count']>r['fresh_derived_direct_sum_floor']
    if 'published_DGLS_shortening_floor' in r:assert r['repository_count']>r['published_DGLS_shortening_floor']
anchors=[{'dimension':d,'research_threshold_to_beat':v,'kind':'published anchor; no registered improvement in RESULTS'} for d,v in [(47,23766960),(48,52416000),(72,6218175600)]]
heads={}
for k in ('packingstar','lin33','lin32','knic','li'):
    a=json.loads(read(HERE/(k+'.txt')))[0];heads[k]={'sha':a['sha'],'commit_utc':a['commit']['committer']['date'],'message':a['commit']['message'].splitlines()[0]}
files=[ROOT/'README.md',ROOT/'RESULTS.md',ROOT/'common/published.py',*sorted((ROOT/'common').glob('source_review*.json'))]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
snapshot={'audit_date_context':'2026-11-04','built_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Fresh selected primary sources plus existing repository. Not an exhaustive literature survey, not a novelty or worldwide priority certification. Counts are existing results, not campaign discoveries.','cohn_scope':'48 rows, dimensions 1–48 only. 25–46 equal repository counts and cite Kravatskiy repository (reference 8); not independent construction verification.','cohn_saved_path':retrieval['saved_path'],'cohn_html_sha256':retrieval['sha256'],'source_heads':heads,'local_source_sha256':manifest,'fresh_cohn_lower':cohn,'anchors':anchors,'active_dimensions5_16':small,'registered_dimensions25_96':rows,'notes':['69 registered rows in 25–96; together with anchors cover all 72 dimensions.','5–16D are equally active targets in the same construction language; only 8D has matching bounds. Existing 25–96D data retained.','Historical comparison columns intentionally retained, not presented as current competitors.','Updated Cohn lower values for 18 and 19 raise derived Gamma72 direct-sum floors in 90 and 91, without threatening repository totals.','No numeric competitor supersession found in selected sources. Above 48 Cohn provides no coverage.','No retrieved commit date is after context audit date; search-engine crawl dates are not publication dates.']}
(HERE/'snapshot.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2)+'\n')
print('PASS: 69 registered rows, 3 anchors; 22 Cohn rows 25–46 equal repo; derived comparisons below records')
print('fragile named margins',[(r['dimension'],r['margin_over_named_competitor']) for r in rows if r.get('margin_over_named_competitor',999999)<=150])
print('fresh Gamma direct sums',[(r['dimension'],r.get('fresh_derived_direct_sum_floor')) for r in rows if r['dimension'] in (90,91)])

print('PASS: all 12 fresh 5–16D lower bounds match repository pins; only8D has matching upper/lower bounds')

print('PASS: all11 saved source bytes match retrieval hashes; Cohn HTML preserved at cohn.html')
