"""Rebuild all weight12 dual supports; audit every blocker-kernel pair.

Uses only the already installed C++ compiler and a bounded temporary binary.
No source data or executable from another task's executor is assumed.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile
from time import monotonic

from audit32_projection import SOURCE,rref,span,require


def main():
    start=monotonic()
    matrix=json.loads(SOURCE.read_text())['generator_matrix']
    G=[sum(x<<j for j,x in enumerate(row)) for row in matrix]
    reduced,pivots=rref(G,32)
    H=[]
    for j in range(32):
        if j in pivots:continue
        H.append((1<<j)|sum(((row>>j)&1)<<p for row,p in zip(reduced,pivots)))
    dual=span(H)
    require(min(h.bit_count() for h in dual if h)==8,'Dual minimum distance')
    supports=sorted(h for h in dual if h.bit_count()==12)
    require(len(pivots)==17 and len(supports)==3584,'Kernel input ranks/counts')
    source=Path(__file__).with_suffix('.cpp')
    with tempfile.TemporaryDirectory(prefix='sol32-kernel-audit-') as temporary:
        executable=Path(temporary)/'audit'
        subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(executable)],check=True)
        payload=f'{len(G)} {len(supports)}\n'+' '.join(map(str,G+supports))+'\n'
        result=subprocess.run([str(executable)],input=payload,text=True,capture_output=True,check=True)
    lines=result.stdout.splitlines()
    histogram={int(k):v for k,v in json.loads(lines[0]).items()}
    pairs=json.loads(lines[1])
    require(histogram=={0:3785280,1:1956448,2:608448,3:69360,4:1200},'Blocker-kernel pair histogram')
    require(sum(histogram.values())==len(supports)*(len(supports)-1)//2==6420736,'Complete pair count')
    require(max(histogram)==4,'Maximum kernel intersection dimension')
    require(len(pairs)==1200,'Four-dimensional intersection pairs')
    kernels={}
    def kernel(index):
        if index not in kernels:
            mask=supports[index]
            constraints=[sum(((g>>j)&1)<<i for i,g in enumerate(G)) for j in range(32) if mask>>j&1]
            reduced,positions=rref(constraints,17)
            basis=[]
            for j in range(17):
                if j not in positions:
                    basis.append((1<<j)|sum(((row>>j)&1)<<p for row,p in zip(reduced,positions)))
            require(len(basis)==6,'Message kernel dimension')
            kernels[index]=basis
        return kernels[index]
    spans=set()
    for i,j in pairs:
        rows,_=rref(kernel(i)+kernel(j),17)
        require(len(rows)==8,'Kernel pair span dimension')
        spans.add(tuple(rows))
    require(len(spans)==1200,'Distinct eight-dimensional pair spans')
    report={'status':'PASS','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'cpp_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'weight12_dual_supports':len(supports),'kernel_dimensions':6,'dual_minimum_distance':8,
            'unordered_kernel_pairs_checked':sum(histogram.values()),
            'kernel_intersection_dimension_histogram':histogram,
            'maximum_affine_blocker_intersection':16,
            'minimum_union_of_two_64_word_blocker_sets':112,
            'distinct_eight_dimensional_pair_spans':len(spans),
            'maximum_kernels_in_any_eight_dimensional_space':2,
            'maximum_minimum_blocker_candidates_freed_by_affine_eight_flat':8,
            'affine_eight_flat_reason':'Three6D kernels in an8D space would give three4D-intersection pairs with the same8D span, contrary to pair-span uniqueness. Each of at most two kernels has four cosets inside an affine8-flat.',
            'reason':'K_S intersection K_T is the kernel of the joint projection on S union T; any nonempty affine intersection has the same direction intersection and at most2^4 words.',
            'scope':'Equal-amplitude norm32 weight12 candidate alphabet over this fixed dense code. Fewer than112 dense deletions can free at most one candidate; candidates with larger blocker sets need at least128. Shared-deletion exchanges of larger size and other alphabets remain open.',
            'new_sota':False,'seconds':monotonic()-start}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
