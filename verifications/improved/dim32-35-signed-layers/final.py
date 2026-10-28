"""Derive all four counts from certificate populations; verify.py proves the angles."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent

def results():
    baselines=json.loads((ROOT/'source/public/baselines.json').read_text())
    for b in baselines:
        file=ROOT/'source/public'/b['source'].rsplit('/',1)[-1]
        assert hashlib.sha256(file.read_bytes()).hexdigest()==b['sha256']
        assert b['kissing_points']==131072+128*b['supports']+2*b['dimension']*(b['dimension']-1)
    out=[]
    for d in range(32,36):
        file=ROOT/('source/multi35/certificate35.json' if d==35 else 'source/classic/certificate%d.json'%d)
        c=json.loads(file.read_text());assert c['dimension']==d and len(c['dense_generator'])==17
        assert c['comparison']==baselines[d-32]
        if d==35:
            assert len(c['removed_root_lines'])==3 and len(c['head_antipodal_lines'])==6
            assert all(type(x) is int and abs(x)<=3 for h in c['head_antipodal_lines'] for x in h)
            total=(1<<17)+128*len(c['full_parity_supports'])+2*d*(d-1)-2*len(c['removed_root_lines'])+2*len(c['head_antipodal_lines'])
        else:
            total=(1<<17)+128*len(c['full_even_supports'])+2*len(c['exceptional_antipodal_lines'])+2*d*(d-1)
            if c.get('heavy_pair_extension'):
                e=c['heavy_pair_extension']
                assert e['removed_root_antipodal_lines']==1 and e['added_head_antipodal_lines']==2
                total+=2*(e['added_head_antipodal_lines']-e['removed_root_antipodal_lines'])
        assert total==c['points']=={32:352800,33:370368,34:390982,35:415314}[d]
        if d in (33,34,35):
            trade=json.loads((ROOT/'root_trades'/('certificate%d.json'%d)).read_text())
            assert trade['format']=='root_star_trade_v1' and trade['dimension']==d
            assert trade['baseline']==file.relative_to(ROOT).as_posix()
            assert hashlib.sha256(file.read_bytes()).hexdigest()==trade['baseline_sha256']
            total+=16*len(trade['stars'])-2*len(trade['removed_root_lines'])
            assert total==trade['points']
        out.append(dict(dimension=d,total=total))
    return out

if __name__=='__main__':
    for r in results():print(json.dumps(r))
