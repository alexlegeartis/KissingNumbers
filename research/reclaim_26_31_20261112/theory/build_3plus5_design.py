"""Reproduce the explicit three-core/five-tail support design."""
import argparse
import itertools
import json
from pathlib import Path

def multiply(a,b):
    result=0
    while b:
        if b&1:result^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return result

def build():
    edges=list(itertools.combinations(range(1,8),2))
    def cube(a):return multiply(multiply(a,a),a)
    first=[a^b for a,b in edges];second=[cube(a)^cube(b)for a,b in edges]
    assert len(set(zip(first,second)))==21
    adjacent=[set()for _ in edges]
    for i,j in itertools.combinations(range(21),2):
        if set(edges[i])&set(edges[j])or first[i]==first[j]or second[i]==second[j]:adjacent[i].add(j);adjacent[j].add(i)
    colors=[-1]*21;nodes=0
    def search(used):
        nonlocal nodes
        nodes+=1
        if min(colors)>=0:return True
        v=max((i for i in range(21)if colors[i]<0),key=lambda i:(len({colors[j]for j in adjacent[i]if colors[j]>=0}),len(adjacent[i])))
        forbidden={colors[j]for j in adjacent[v]if colors[j]>=0}
        for color in range(min(used+1,8)):
            if color not in forbidden:
                colors[v]=color
                if search(max(used,color+1)):return True
        colors[v]=-1;return False
    assert search(0)
    triples=[[first[i]-1,8+second[i]-1,16+colors[i]]for i in range(21)]
    for i,j in itertools.combinations(range(21),2):assert len(set(triples[i])&set(triples[j]))<=1-int(bool(set(edges[i])&set(edges[j])))
    return {'description':'Three core coordinates per five-coordinate tail support, with all 21 omitted edges of K7. First two colors use GF(8) APN x^3; third is exact graph coloring.','tail_omitted_edges':edges,'core_triples':triples,'colors_used':len(set(colors)),'graph_search_nodes':nodes}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args();text=json.dumps(build(),indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')
