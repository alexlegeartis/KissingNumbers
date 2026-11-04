"""Independent standard-library audit of the supplied 32D projection claim.

This certifies facts about one fixed dense code, not a kissing-number upper
bound or an obstruction to coupled deletions and insertions.
"""
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'verifications/improved/dim32-35-signed-layers/source/public/dense_code.json'


def require(value,message):
    if not value:raise ValueError(message)


def rref(rows,n):
    values=list(rows);pivots=[];rank=0
    for column in range(n):
        index=next((i for i in range(rank,len(values)) if values[i]>>column&1),None)
        if index is None:continue
        values[rank],values[index]=values[index],values[rank]
        for i in range(len(values)):
            if i!=rank and values[i]>>column&1:values[i]^=values[rank]
        pivots.append(column);rank+=1
    return values[:rank],pivots


def span(rows):
    result=[0]
    for row in rows:result += [x^row for x in result]
    require(len(set(result))==len(result),'Dependent span basis')
    return result


def main():
    data=json.loads(SOURCE.read_text());matrix=data['generator_matrix']
    require(len(matrix)==17 and all(len(row)==32 and all(type(x) is int and x in (0,1) for x in row) for row in matrix),'Binary source shape')
    G=[sum(x<<j for j,x in enumerate(row)) for row in matrix]
    reduced,pivots=rref(G,32);require(len(pivots)==17,'Parent rank')
    free=[j for j in range(32) if j not in pivots]
    H=[]
    for j in free:
        row=1<<j
        for value,pivot in zip(reduced,pivots):
            if value>>j&1:row|=1<<pivot
        H.append(row)
    require(len(H)==15 and all((h&g).bit_count()%2==0 for h in H for g in G),'Dual rank and orthogonality')
    dual=span(H);words=span(G)
    require(min(w.bit_count() for w in words if w)==8,'Parent distance')
    histogram=dict(sorted(Counter(w.bit_count() for w in dual).items()))
    expected={0:1,8:124,10:1152,12:3584,14:6016,16:11014,18:6016,20:3584,22:1152,24:124,32:1}
    require(histogram==expected,'Dual histogram differs')
    weight8=[w for w in dual if w.bit_count()==8]
    rank10={a|b for a,b in combinations(weight8,2) if (a|b).bit_count()==12}
    require(len(rank10)==264,'Rank-ten support count')
    for mask in rank10:
        inside=[h for h in weight8 if h&mask==h]
        require(len(inside)==3 and inside[0]^inside[1]==inside[2],'Rank-ten shortening')
    counts={'rank10':len(rank10),'rank11_weight8':124*comb(24,4)-3*len(rank10),
            'rank11_weight10':1152*comb(22,2),'rank11_weight12':3584}
    counts['rank12']=comb(32,12)-sum(counts.values())
    require(counts=={'rank10':264,'rank11_weight8':1316832,'rank11_weight10':266112,'rank11_weight12':3584,'rank12':224206048},'Support census')
    # Plotkin's elementary weight sum: (2^k-1)*8 <= 2^(k-1)*s.
    require((2**4-1)*8>2**3*14,'Dimension-four exclusion through support size14')
    require((2**3-1)*8>2**2*12,'Dimension-three exclusion through support size12')
    # At k=3,s=14 equality forces all seven nonzero words to weight8.
    # The following seven column types, each occurring twice, solve the seven
    # weight equations. Character orthogonality gives uniqueness of this
    # solution (the inverse 8-point Walsh transform).
    require(all(sum(2*((u&v).bit_count()%2) for v in range(1,8))==8 for u in range(1,8)),'Dimension-three equality profile')
    # Select actual supports of each shortening type, then enumerate all
    # parent words and all4096 sign patterns independently of the formulas.
    small=[h for h in dual if 0<h.bit_count()<=12]
    def shortening(mask):return [h for h in small if h&mask==h]
    examples={'rank10':min(rank10)}
    for weight in (8,10,12):
        found=None
        for h in (h for h in small if h.bit_count()==weight):
            outside=[j for j in range(32) if not h>>j&1]
            for extra in combinations(outside,12-weight):
                mask=h|sum(1<<j for j in extra)
                if len(shortening(mask))==1:found=mask;break
            if found is not None:break
        require(found is not None,'Missing rank-eleven example')
        examples[f'rank11_weight{weight}']=found
    for indices in combinations(range(32),12):
        mask=sum(1<<j for j in indices)
        if not shortening(mask):examples['rank12']=mask;break
    sample_reports={}
    expected_blockers={'rank10':{128:1024,512:3072},'rank11_weight8':{320:2048,512:2048},
                       'rank11_weight10':{192:2048,640:2048},'rank11_weight12':{64:2048,768:2048},
                       'rank12':{416:4096}}
    for label,mask in examples.items():
        positions=[j for j in range(32) if mask>>j&1]
        frequencies=Counter(word&mask for word in words)
        distribution=Counter()
        for signs in range(1<<12):
            pattern=sum(((signs>>i)&1)<<j for i,j in enumerate(positions))
            blockers=frequencies[pattern]+sum(frequencies[pattern^(1<<j)] for j in positions)
            distribution[blockers]+=1
        require(dict(distribution)==expected_blockers[label],'Explicit parent enumeration differs from algebraic blocker counts')
        sample_reports[label]={'support_mask':mask,'shortened_dual_nonzero_words':shortening(mask),
                               'projection_size':len(frequencies),'blocker_histogram':dict(distribution)}
    # Equal-amplitude norm32 support-s vectors: one mismatch still violates
    # dot<=16 iff (s-2)^2>8*s. No approximate square roots decide this.
    require(all((s-2)**2>8*s for s in (12,13,14)),'Sparse normalized conflict threshold')
    certificate=ROOT/'verifications/improved/dim32-35-signed-layers/source/classic/certificate32.json'
    sparse=json.loads(certificate.read_text())
    require(rref(sparse['dense_generator'],32)[0]==reduced,'Certificate dense code is the same coordinate code')
    octads=sparse['full_even_supports']
    require(len(octads)==1704 and len(set(octads))==1704 and all(o.bit_count()==8 for o in octads),'Complete even-sign octad supports')
    # For overlap t, dot=2*sqrt(8/3)*(t-2r)>16 iff (t-2r)^2>24
    # and t-2r>0. Enumerate all128 allowed even-sign patterns exactly, for
    # both possible target sign parities; outside-intersection signs are free.
    per_overlap={}
    for t in range(9):
        values=[]
        for target in (0,1):
            count=0
            for signs in range(256):
                if signs.bit_count()%2:continue
                disagreement=((signs^target)&((1<<t)-1)).bit_count()
                z=t-2*disagreement
                if z>0 and z*z>24:count+=1
            values.append(count)
        per_overlap[t]=min(values)
    require(per_overlap=={0:0,1:0,2:0,3:0,4:0,5:4,6:2,7:8,8:1},'Even-sign blocker multiplicities')
    overlap_counts=[];sparse_prices=[]
    for support in (h for h in dual if h.bit_count()==12):
        histogram_s=Counter((support&o).bit_count() for o in octads)
        overlap_counts.append(sum(v for t,v in histogram_s.items() if t>=5))
        sparse_prices.append(sum(v*per_overlap[t] for t,v in histogram_s.items()))
    require(min(overlap_counts)==140 and min(sparse_prices)==498,'Minimum octad blocker prices')
    report={'status':'PASS','source':str(SOURCE.relative_to(ROOT)),
            'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'parent_rank':17,'parent_distance':8,'dual_rank':15,'dual_distance':8,
            'dual_weight_histogram':histogram,'weight12_support_census':counts,
            'explicit_blocker_examples':sample_reports,
            'projection_radius_through_size14':'at most1, by dual shortening/Plotkin and nonzero syndrome-column coverage',
            'minimum_dense_blockers_equal_amplitude_support12':64,
            'complete_octad_certificate':str(certificate.relative_to(ROOT)),
            'complete_octad_certificate_sha256':hashlib.sha256(certificate.read_bytes()).hexdigest(),
            'full_even_octad_supports':len(octads),
            'minimum_full_octad_supports_with_overlap_at_least5':min(overlap_counts),
            'minimum_even_signed_octad_blockers':min(sparse_prices),
            'minimum_total_dense_plus_complete_octad_blockers':64+min(sparse_prices),
            'minimum_octad_blockers_by_overlap':per_overlap,
            'prior_history':'KNOWLEDGE section130 already discusses a16-coordinate radius2 projection obstruction; section132 has exact dense covering radius6,35-optimality and single-coset injection. This refines fixed-layer obstructions, not a new construction.',
            'scope':'Fixed full dense layer only. Does not exclude shared-deletion exchanges, replacing the dense code or any kissing-number improvement.',
            'new_sota':False}
    target=Path(__file__).with_name('audit32_projection.json');target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','dual_distance':8,'rank10_supports':264,'minimum_weight12_dense_blockers':64}))


if __name__=='__main__':main()
