"""Compare near-minimum cap counts with literal pairs of partial Dyck words.

Boundary order is outer-to-inner on the upper side, then inner-to-outer below.
In the literal calculation every physical vertex belongs to an explicit arc-graph
component. Graphs containing a closed component contribute zero.
"""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import cache
from hashlib import sha256
from pathlib import Path
import json

from minimum_caps import matchings, partial_words


def literal_components(top,lower):
    upper_arcs,T=top
    lower_arcs,B=lower
    n=len(upper_arcs)
    parent=list(range(n))

    def find(i):
        while i!=parent[i]:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i

    for arcs in (upper_arcs,lower_arcs):
        for i,j in enumerate(arcs):
            if j>i:
                parent[find(i)]=find(j)
    sizes=Counter(find(i) for i in range(n))
    bound=defaultdict(list)
    for index,vertex in enumerate((*T,*reversed(B))):
        bound[find(vertex)].append(index)
    if len(bound)!=len(sizes):
        return None
    mates=[-1]*(len(T)+len(B))
    exceptional=0
    for root,ends in bound.items():
        assert len(ends)==2
        i,j=ends
        mates[i]=j
        mates[j]=i
        shortest=1 if i<len(T)<=j else 2
        extra=sizes[root]-shortest
        assert extra>=0 and extra%2==0
        exceptional+=extra>0
    return tuple(mates),exceptional


@cache
def literal_histogram(a,b,t):
    upper=[u for u in partial_words(t) if len(u[1])==a]
    lower=[v for v in partial_words(t) if len(v[1])==b]
    hist=Counter()
    for u in upper:
        for v in lower:
            result=literal_components(u,v)
            if result is not None:
                hist[result]+=1
    return hist,len(upper)*len(lower)


def schedule():
    cases=[]
    for k in range(1,5):
        for p in matchings(k):
            for a in range(2*k+1):
                j=sum(x>=a for x in p[:a])
                L=2*k-j
                for d,bound in ((1,7),(2,6)):
                    if L+2*d<=bound:
                        cases.append((a,p,d))
    # (13,0),(1,2),...,(11,12) has two through paths and five upper pairs.
    # Target n here denotes the number of exposed ends, 2k; t=L+2=14.
    ends=14
    shifted=tuple((i+1)%ends if i%2 else (i-1)%ends for i in range(ends))
    cases.append((12,shifted,1))
    return cases


def run(candidate=None):
    cases=schedule()
    assert len(cases)==len(set(cases))
    records=[];small_by_defect=Counter();nonempty_all_exceptional=Counter()
    ranges={}
    for a,p,d in cases:
        b=len(p)-a
        through=sum(x>=a for x in p[:a])
        L=len(p)-through
        t=L+2*d
        hist,pair_count=literal_histogram(a,b,t)
        ranges[a,b,t]=pair_count
        distribution=[hist[p,j] for j in range(d+1)]
        value=sum(distribution)
        assert hist[p,0]==0
        if candidate is not None:
            result=candidate(a,p,d)
            assert type(result) is int and result==value,(a,p,d,value,result)
        if len(p)<=8:
            small_by_defect[d]+=1
            if len(p)//2<=d:
                nonempty_all_exceptional[d]+=distribution[len(p)//2]>0
        else:
            assert (a,b,through,L,t)==(12,2,2,12,14)
            assert value==(len(p)//2)**2==49
        records.append(dict(a=a,b=b,p=list(p),through=through,L=L,d=d,t=t,
                            distribution_by_number_of_exceptional_components=distribution,
                            count=value))
    assert all(nonempty_all_exceptional[d]>0 for d in (1,2))
    assert all(any(r["a"]==a and r["b"]==b and r["d"]==1 for r in records)
               for a,b in ((0,2),(2,0),(1,1),(2,2)))
    return dict(smaller_case_counts={str(d):small_by_defect[d] for d in (1,2)},
                small_types_admitting_no_ordinary_component={str(d):nonempty_all_exceptional[d] for d in (1,2)},
                literal_word_pairs=sum(ranges.values()),literal_height_time_groups=len(ranges),
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),checks=records)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--write",type=Path)
    parser.add_argument("--literal-only",action="store_true")
    options=parser.parse_args()
    if options.literal_only:
        if options.write:parser.error("A saved cross-check requires the candidate implementation")
        data=run()
    else:
        from excess_caps import count_near_minimum
        data=run(count_near_minimum)
        code=Path(__file__).resolve().with_name("excess_caps.py")
        data["candidate_sha256"]=sha256(code.read_bytes()).hexdigest()
        if options.write:
            tmp=options.write.with_suffix(".pending")
            tmp.write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")
            tmp.replace(options.write)
        else:
            local=code.parent/"certificates/excess-caps.json"
            if not local.is_file():
                local=code.parents[3]/"output/research73/star/excess-caps.json"
            assert json.loads(local.read_bytes())==data
    print(data["smaller_case_counts"],data["small_types_admitting_no_ordinary_component"],
          data["literal_word_pairs"],data["checks"][-1]["count"])
