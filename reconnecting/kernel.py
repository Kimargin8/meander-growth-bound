"""Exact, independently bounded first-return weights in a bounded physical strip."""
from fractions import Fraction as Q
from functools import cache
from math import prod

DEN = 4096
DEPTH = 6
DEGREE = 83
TARGET = Q(12)


def down(x,den):
    return Q(x.numerator*den//x.denominator,den)


def path_returns(width, edges, x, maximum=3):
    # D_j = det(I - A_j/s). Interior embeddings are symmetric under reflection.
    assert tuple(edges) == tuple(reversed(edges)) and len(edges)==width+2
    det = [Q(1),Q(1)]
    for edge in edges:
        det.append(det[-1] - edge*edge*det[-2]/x)
    assert all(d>0 for d in det)
    size = width+3
    by_delta = {}
    for delta in range(1,maximum+1):
        row = []
        for u in range(width-delta+1):
            a = u+1
            b = a+delta
            row.append(prod(edges[a:b])*det[a]*det[size-1-b]/det[size]/
                       x**(1+delta//2))
        by_delta[delta] = row
    return by_delta


def continuants(width,depth,x):
    finer = 2**30
    deeper = [Q(1)]*(width+2*depth)
    for level in reversed(range(depth)):
        values=path_returns(width+2*level,deeper,x,3 if level==0 else 1)
        if level:
            deeper=[1+down(k,finer) for k in values[1]]
    return values


@cache
def literal(width, depth, degree):
    """Closed coefficient arrays indexed by physical (single-site) word length."""
    child = [[0,1]+[0]*(degree-1) for _ in range(width+2)]
    if depth>1:
        inside = literal(width+2, depth-1, degree-2)[1]
        child = [list(row) for row in inside]
        for row in child: row[1] = 1
    out = {d:[[0]*(degree+1) for _ in range(width+1-d)] for d in (1,2,3)}
    for u in range(width):
        a=u+1
        counts=[[0]*(width+3) for _ in range(degree-1)]
        counts[0][a]=1
        for length, layer in enumerate(counts):
            for d,rows in out.items():
                if u < len(rows):
                    rows[u][length+2]=layer[a+d]
            for p,count in enumerate(layer):
                if count:
                    for edge,step in ((p,1),(p-1,-1)):
                        if 0<=edge<len(child):
                            seq=child[edge]
                            for delta in range(1,min(degree-1-length,len(seq)),2):
                                if seq[delta]:
                                    counts[length+delta][p+step] += count*seq[delta]
    return out


def tables(w,x=TARGET):
    poly = literal(w,DEPTH,DEGREE)
    species = continuants(w,DEPTH,x)
    a = {d:[int(DEN*t) for t in row] for d,row in species.items()}
    b = {}
    for d,rows in poly.items():
        b[d]=[int(DEN*sum(Q(v,x**(i//2)) for i,v in enumerate(seq) if v)) for seq in rows]
    return a,b
