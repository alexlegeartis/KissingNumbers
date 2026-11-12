"""Fresh exact replay of all sign classes, support interactions and transfers."""
import hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from verify_family import check
from geometry_counts import verify_geometry,counts
from section63 import verify_section
from expanded_geometry import check_expanded
from replay_parents import verify_parents
from final import results

def integrity():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    manifest=json.loads((HERE/'SHA256.json').read_text())
    actual={p.relative_to(HERE).as_posix() for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='SHA256.json'}
    assert set(manifest)==actual,'manifest coverage'
    for name,digest in manifest.items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    inputs=json.loads((HERE/'data/source_inputs.json').read_text())
    repository=HERE.parents[2]
    for name,digest in inputs['sha256'].items():assert hashlib.sha256((repository/name).read_bytes()).hexdigest()==digest,name
    return dict(files=len(manifest),pinned_source_files=len(inputs['sha256']))

def verify():
    checked=integrity()
    selection=json.loads((HERE/'data/selection.json').read_text())
    names=sorted({selection[k] for k in ('short','global','expanded66')})
    families={name:check(HERE/name,cpp=True) for name in names}
    from verify_independent import check as independent_check
    independent=independent_check()
    geometry=verify_geometry();parents=verify_parents()
    section=verify_section(HERE/selection['global'],enumerate_words=True)
    assert independent['section63']['total']==section['total']
    expanded=check_expanded(HERE/selection['expanded66'],66,cpp=True)
    assert expanded['terms']['weight4']==parents['four_support_counts'][66]*16
    rows=results();totals={r['dimension']:r['total'] for r in rows}
    assert totals[63]==section['total'] and totals[66]==expanded['total']
    assert all(totals[int(d)]==value for d,value in independent['totals'].items())
    assert all(r['gain']>0 for r in rows)
    subcode=None
    if 'subcode67' in selection:
        from baseline67 import verify_baseline67
        from subcode67 import verify67
        from verify_subcode67_independent import review67
        baseline=verify_baseline67()
        subcode=verify67(HERE/selection['subcode67'])
        direct=review67(HERE/selection['subcode67'])
        assert direct['total']==subcode['total'] and direct['certificate_sha256']==subcode['certificate_sha256']
        subcode['independent_direct_review']={k:v for k,v in direct.items() if k!='high_overlap_pairs'}
        assert baseline['expanded']['terms']['weight4']==parents['four_support_counts'][67]*16
        subcode['baseline_replay']=baseline
        assert totals[67]==subcode['total']
    return dict(status='PASS',integrity=checked,dimensions=rows,families=families,
                independent_full_family_review=independent,
                geometry=geometry,parent_replays=parents,section63=section,expanded66=expanded,subcode67=subcode)

if __name__=='__main__':print(json.dumps(verify(),indent=2))
