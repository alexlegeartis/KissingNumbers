from collections import Counter
def span(generators):
    h = {0}
    for g in generators:
        h |= {w ^ g for w in h}
    return h
p = sum(1 << j for j in (0,1,5,6,7,9,11))
g = [(p << j) | (1 << 23) for j in range(12)]
words = {w & ((1 << 21)-1) for w in span(g)}
lines = sorted(w for w in words if w.bit_count() == 5)
assert len(words) == 4096 and len(lines) == 21 and len(span(lines)) == 1024
coord = {0: 0}
for l in lines:
    z = l ^ lines[0]
    if z not in coord:
        b = len(coord)
        coord.update({w ^ z: i ^ b for w, i in list(coord.items())})
assert len(coord) == 512
translated = [coord[l ^ lines[0]] for l in lines]
weights = Counter(sum((a & t).bit_count() % 2 for t in translated) for a in range(512))
assert weights == {0:1, 5:16, 8:130, 9:160, 12:120, 13:80, 16:5}
both = set(weights) | {21-w for w in weights}
assert set(range(8,14)) & both == {8,9,12,13}

C = span(lines)
masks = (317651,405196,1307954,498987)
K0 = {w for w in C if all((w&m).bit_count()%2 == 0 for m in masks)}
U = span((1185854,1258025,176283,236172,1005402,949581,2020840))
v = 136264
assert len(K0)==64 and len(U)==16 and U<=K0 and v not in K0
H = (K0-U) | {v^w for w in U}
assert len(H)==64 and all(w^p not in H for w in H for p in lines)
for b in (0,1):
    h = [w for w in H if w.bit_count()%2==b]
    assert len(h)==32 and len(span([w^h[0] for w in h]))==64
hist = Counter(sum(w^p in H for p in lines) for w in C-H)
assert hist == {1:768,3:192}
lab = {w:next(p for p in lines if w^p in H) for w in C-H}
assert all(len({lab[w],lab[w^p],p})==3 for w in lab for p in lines if w^p in lab)
