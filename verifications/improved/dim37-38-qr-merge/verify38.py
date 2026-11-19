#!/usr/bin/env python3
"""Exact replay of the 38-dimensional merged-QR layered configuration."""
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

P = Path(__file__).resolve().parent
SOURCE_SHA = "ccb530f90e914df84cbc8d19a71d30c79146d23a2e8e86e7a813d4d92a5a6a8b"

def projected_dense(word):
    a, b = word & 1, word >> 1 & 1
    if (a, b) == (1, 0):
        return None
    return (3 * (1 - 2*a) if a == b else 0, *(
        2 * (1 - 2 * ((word >> j) & 1)) for j in range(2, 39)
    ))

def no_collision(d, n1, n2):
    return d <= 0 or 4*d*d <= n1*n2

def main():
    data = json.loads((P / "dense_parent39.json").read_text())
    G = [int(g, 16) for g in data["generators_hex"]]
    R = [int(r, 16) for r in data["representatives_hex"]]
    assert len(G) == 15 and len(R) == len(set(R)) == 10
    assert all(0 <= w < (1<<39) for w in G+R)
    H = [0]
    for g in G:
        H += [h^g for h in H]
    assert len(H) == len(set(H)) == 32768
    min_kernel = min(h.bit_count() for h in H if h)
    assert min_kernel == 12
    distances = [min((x^y^h).bit_count() for h in H)
                 for x,y in itertools.combinations(R,2)]
    assert distances == [10]*45
    parent = Counter(x & 3 for r in R for h in H for x in [r^h])
    assert parent == {0:81920, 1:81920, 2:81920, 3:81920}
    # Low-2-bit 0/3 supply the equal class, 2 supplies the fixed opposite
    # orientation (+,-). Each class is injective: a nonzero merged entry
    # determines both old signs; the zero entry determines the fixed pair.
    E, O = parent[0]+parent[3], parent[2]
    assert E == 163840 and O == 81920
    assert projected_dense(0) == (3,)+(2,)*37
    assert projected_dense(1) is None
    assert projected_dense(2) == (0,)+(2,)*37
    assert projected_dense(3) == (-3,)+(2,)*37

    # All distinct parent products are at most 39−2*10=19.
    # After multiplying by four and merging, the product increment is
    # ≤1 (equal/equal), 0 (mixed), or −8 (both fixed opposite).
    # The sparse potentials only use maximum coordinate magnitudes.
    maxima = {
        "equal_equal": (77,157,157), "equal_opposite": (76,157,148),
        "opposite_opposite": (68,148,148), "equal_octad": (34,157,32),
        "opposite_octad": (32,148,32), "equal_root": (20,157,32),
        "opposite_root": (16,148,32), "octad_octad": (16,32,32),
        "octad_root": (16,32,32), "root_root": (16,32,32),
    }
    assert all(no_collision(d,n1,n2) for d,n1,n2 in maxima.values())

    source = (P / "supports_3088.txt").read_bytes()
    assert hashlib.sha256(source).hexdigest() == SOURCE_SHA
    rows = source.decode("ascii").splitlines()
    assert len(rows) == 3088
    assert all(len(s) == 38 and set(s) <= {"0","1"} and s.count("1") == 8 for s in rows)
    S = [int(s[::-1],2) for s in rows]
    assert len(set(S)) == 3088
    intersections = Counter((s&t).bit_count() for s,t in itertools.combinations(S,2))
    assert intersections == {0:562756,1:1540993,2:1682649,3:710496,4:269434}
    # Same octad: even parity ensures sign distance >=2. Distinct
    # octads meet in <=4. Distinct roots meet in one place or have
    # a sign difference on their common two coordinates.
    M, Q = 128*len(S), 2*38*37
    assert M == 395264 and Q == 2812
    total = E+O+M+Q
    assert total == 643836
    report = {
        "status":"PASS", "dimension":38, "points":total,
        "dense_equal":E, "dense_fixed_opposite":O, "middle":M, "roots":Q,
        "min_parent_kernel_distance":min_kernel,
        "cross_coset_pairs_at_distance_10":len(distances),
        "merged_parent_coordinates_zero_based":[0,1],
        "middle_support_sha256":SOURCE_SHA,
        "octad_intersections":dict(sorted(intersections.items())),
        "integer_product_ceiling_norms":maxima,
    }
    print(json.dumps(report, indent=2))
    return report

if __name__ == "__main__":
    main()
