"""Reproduce the certificate exactly from the retained witness and explicit trades."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    source=HERE.parent/'phase/certificate.json'
    delta=json.loads((HERE/'delta.json').read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest()==delta['baseline_sha256']
    c=json.loads(source.read_text())
    old_phase=set(map(tuple,c['phase_orbits']))
    removed=set(map(tuple,delta['phase_removed']));added=set(map(tuple,delta['phase_added']))
    assert removed<=old_phase and not added.intersection(old_phase)
    new_phase=(old_phase-removed)|added
    c['phase_orbits']=[list(row) for row in sorted(new_phase)]
    roots=delta['root_labels']
    assert len(roots)==len(set(roots))==252 and roots==sorted(roots)
    index={label:j for j,label in enumerate(roots)}
    old_lifts=set()
    for cell in c.pop('lift_body_cosets'):
        nonzero=[j for j,x in enumerate(cell) if x]
        for p in range(3):
            labels=[index[x] if x else 252 for x in cell]
            phases=[0,0,0];phases[nonzero[1]]=p
            old_lifts.add(tuple(labels+phases[1:]))
    removed=set(map(tuple,delta['lift_removed']));added=set(map(tuple,delta['lift_added']))
    assert removed<=old_lifts and not added.intersection(old_lifts)
    sorted_lifts=sorted((old_lifts-removed)|added)
    c['lift_phase_orbits']=[list(row) for row in sorted_lifts]
    c['format']='k12_phase_and_lift_orbits_v2'
    c['points_36']=403704+3*len(c['phase_orbits'])
    c['points_37']=c['points_36']+6*len(c['lift_phase_orbits'])
    # JSON object order is immaterial mathematically; match the shipped serialization.
    target=json.loads((HERE/'certificate.json').read_text())
    assert c==target
    raw=(json.dumps({key:c[key] for key in target},separators=(',',':'))+'\n').encode()
    assert raw==(HERE/'certificate.json').read_bytes()
    print(json.dumps(dict(reproduced=True,points_36=c['points_36'],points_37=c['points_37'],
                          certificate_sha256=hashlib.sha256(raw).hexdigest())))


if __name__=='__main__':main()
