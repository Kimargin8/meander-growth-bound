from fractions import Fraction as Q
from kernel import path_returns
DEN=4096
X=Q(1173,100)
DEPTH=8
ROUNDS=40
PREC=2**36
TERMS=52


def ceil(a,b):
    return -(-a//b)


def ceil_ratio(x,scale):
    return ceil(x.numerator*scale,x.denominator)


def upper_walks(width,edges):
    # edges are integer multiples of 1/PREC. Advance a local mass distribution
    # by two steps and divide by x. Round every entry upward after each step.
    m=width+3
    upper_edge=Q(max(edges),PREC)
    q=4*upper_edge**2/X
    rest=2*upper_edge/X*q**TERMS/(1-q)
    assert q<1 and upper_edge<Q(117,100)
    def advance(vec,numer,denom):
        return [ceil(numer*((edges[j-1]*vec[j-1] if j else 0)+
                            (edges[j]*vec[j+1] if j<m-1 else 0)),denom*PREC) for j in range(m)]
    out=[]
    for u in range(width):
        v=[0]*m;v[u+1]=PREC
        ans=0
        for _ in range(TERMS):
            first=advance(v,X.denominator,X.numerator)
            ans+=first[u+2]
            v=advance(first,1,1)
        out.append(Q(ans,PREC)+rest)
    return out


def weights(method):
    unit=ceil_ratio(Q(7,6),PREC)
    edges=[unit]*(32+2*DEPTH)
    for lev in reversed(range(DEPTH)):
        w=32+2*lev
        if method=="walks":
            bonus=upper_walks(w,edges)
        else:
            bonus=path_returns(w,[Q(v,PREC) for v in edges],X,1)[1]
        if lev:
            edges=[PREC+ceil_ratio(q,PREC) for q in bonus]
    return [ceil_ratio(q,DEN) for q in bonus]

