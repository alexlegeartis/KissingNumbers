from itertools import combinations,product
def hanani(n):
    assert n>1 and n&(n-1)==0
    k=n-1
    B=[t for t in combinations(range(n),4) if t[0]^t[1]^t[2]^t[3]==0]
    def p(b,x): return (b%3)*k+x
    X=[]
    for t in B:
        if k not in t:
            for a,b,c in product(range(3),repeat=3):
                X.append(tuple(p(j,x) for j,x in zip((a,b,c,-a-b-c),t)))
        else:
            t=t[:3]
            for i in range(3):
                for w in t:
                    u,v=[x for x in t if x!=w]
                    X.append((p(i,u),p(i,v),p(i+1,w),p(i+2,w)))
            for a,b in product(range(3),repeat=2):
                X.append((3*k,p(a,t[0]),p(b,t[1]),p(-a-b,t[2])))
    X += [(p(i,u),p(i,v),p(i+1,u),p(i+1,v)) for i in range(3) for u,v in combinations(range(k),2)]
    X += [(3*k,p(0,u),p(1,u),p(2,u)) for u in range(k)]
    X = sorted(tuple(sorted(row)) for row in X)
    assert all(len(set(row))==4 for row in X)
    triples=[t for row in X for t in combinations(row,3)]
    assert len(triples)==len(set(triples)) and set(triples)==set(combinations(range(3*k+1),3))
    return X
if __name__=='__main__':
    X=[row for row in hanani(16) if max(row)<42]
    T=set(t for x in X for t in combinations(x,3))
    C=[x for x in combinations(range(42),4) if all(t not in T for t in combinations(x,3))]
    D=[x for x in combinations(range(42),5) if all(t not in T for t in combinations(x,3))]
    print(len(X),len(C),len(D),C[:20],D[:20])
    for x in C:
       if any(t in T for t in combinations(x,3)):continue
       X.append(x)
       T.update(combinations(x,3))
    print(len(X))

def double(X,n):
    out=[tuple(2*i+v for i,v in zip(q,b)) for q in X for b in product(range(2),repeat=4) if sum(b)%2==0]
    out += [(2*i,2*i+1,2*j,2*j+1) for i,j in combinations(range(n),2)]
    out.sort()
    assert all(len(set(x))==4 for x in out)
    triples=[t for x in out for t in combinations(x,3)]
    assert len(triples)==len(set(triples)) and set(triples)==set(combinations(range(2*n),3))
    return out
