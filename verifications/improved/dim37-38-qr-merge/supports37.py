"""Regenerate Brouwer's 2,742 weight-eight supports on 36 ordinary sites.

Numerical source: https://aeb.win.tue.nl/codes/cwc/d8/a36.8.8.2742H
The 2,739-word part is eight orbits of the diagonal PSL(2,11) action
on three projective lines. Three further masks are appended.
"""
from collections import deque
from hashlib import sha256
SOURCE = "https://aeb.win.tue.nl/codes/cwc/d8/a36.8.8.2742H"
SEEDS = (
    (0x0000035B1,132), (0x000027E10,165), (0x003024115,660),
    (0x003115404,660), (0x00F081104,660), (0x0270001C8,165),
    (0x0271C8000,165), (0x07D101000,132),
)
EXTRAS = (0x0000FF000,0x000FF0000,0x000F0F000)
ORDER_SHA = "75de16ebd0f473d39d851a6daf367b3e0ea00478984513f9ccc68ef7de5405fa"

def supports():
    # Projective labels 0,1,...,10,∞ are stored in positions 0,...,9,11,10.
    r = [*range(10),11,10]
    g = [[0]*12 for _ in range(2)]
    for j,i in enumerate(r):
        g[0][i]=r[(j+1)%11] if j<11 else r[11]
        g[1][i]=r[11] if j==0 else r[0] if j==11 else r[-pow(j,-1,11)%11]
    perm = [tuple(1<<(i+12*k) for k in range(3) for i in p) for p in g]
    all_rows=set()
    for seed,expected in SEEDS:
        orbit={seed}
        Q=deque([seed])
        while Q:
            x=Q.popleft()
            bits=[i for i in range(36) if x>>i&1]
            for p in perm:
                y=sum(p[i] for i in bits)
                if y not in orbit:
                    orbit.add(y);Q.append(y)
        assert len(orbit)==expected and all_rows.isdisjoint(orbit)
        all_rows|=orbit
    assert len(all_rows)==2739 and all_rows.isdisjoint(EXTRAS)
    all_rows.update(EXTRAS)
    result=sorted(all_rows)
    assert sha256("".join(f"{x:09X}\n" for x in result).encode()).hexdigest()==ORDER_SHA
    return result

if __name__=='__main__':
    for row in supports():
        print(f"{row:09X}")
