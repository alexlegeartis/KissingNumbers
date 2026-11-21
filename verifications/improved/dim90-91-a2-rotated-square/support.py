"""Hanani 3n−2 expansion with n=16 and distinguished point 15."""
from itertools import combinations,product
def supports():
    B=[t for t in combinations(range(16),4) if t[0]^t[1]^t[2]^t[3]==0]
    def p(b,x): return (b%3)*15+x
    X=[]
    for t in B:
        if 15 not in t:
            for a,b,c in product(range(3),repeat=3):
                X.append(tuple(p(k,x) for k,x in zip((a,b,c,-a-b-c),t)))
        else:
            t=t[:3]
            for i in range(3):
                for w in t:
                    u,v=[x for x in t if x!=w]
                    X.append((p(i,u),p(i,v),p(i+1,w),p(i+2,w)))
            for a,b in product(range(3),repeat=2):
                X.append((45,p(a,t[0]),p(b,t[1]),p(-a-b,t[2])))
    X += [(p(i,u),p(i,v),p(i+1,u),p(i+1,v)) for i in range(3) for u,v in combinations(range(15),2)]
    X += [(45,p(0,u),p(1,u),p(2,u)) for u in range(15)]
    assert len(X)==3795 and all(len(set(row))==4 for row in X)
    triples=[z for row in X for z in combinations(sorted(row),3)]
    assert len(set(triples))==len(triples)==15180 and set(triples)==set(combinations(range(46),3))
    return [tuple(sorted(row)) for row in X if 45 not in row]
if __name__=='__main__':
    print(len(supports()))
