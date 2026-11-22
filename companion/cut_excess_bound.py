"""Finite upper bounds for a near-shortest side of a connected closed meander."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from math import comb
from pathlib import Path
import json

from check_excess_caps import literal_components
from height_envelope import automaton, forbidden_words, internal_cycle
from minimum_caps import gf_layers, partial_words


def single_factor_opposite(limit):
    """Ballot paths on both banks, excluding the isolated consecutive two-point cycle."""
    layer = {(0, 0, False): 1}
    totals = [Counter({(0, 0): 1})]
    for t in range(1, limit+1):
        nxt = Counter()
        for (a, b, double_open), value in layer.items():
            for x, y in ((1,1), (1,-1), (-1,1), (-1,-1)):
                if a+x < 0 or b+y < 0 or (double_open and x == y == -1):
                    continue
                nxt[a+x, b+y, x == y == 1] += value
        counts = Counter()
        for (a, b, _), value in nxt.items():
            counts[a,b] += value
        totals.append(counts)
        layer = nxt
    return totals


def possible_opposite(limit,rows):
    """Ballot words ending at every ordered height pair in the avoidance DFA."""
    layer={(0,0,0):1}
    totals=[Counter({(0,0):1})]
    for t in range(1,limit+1):
        nxt={}
        for (a,b,state),weight in layer.items():
            for letter,dest in enumerate(rows[state]):
                if dest<0:continue
                h=a+(1 if not letter&2 else -1)
                k=b+(1 if not letter&1 else -1)
                if h<0 or k<0:continue
                key=(h,k,dest)
                nxt[key]=nxt.get(key,0)+weight
        counts=Counter()
        for (a,b,_),weight in nxt.items():counts[a,b]+=weight
        assert all(counts[a,b]==counts[b,a] for a,b in counts)
        totals.append(counts)
        layer=nxt
    return totals


def exceptional_shapes(D):
    # Retain exactly the halves in which every component exceeds its own minimum.
    result = Counter({(0,0,0,0): 1})
    for e in range(1,4*D+1):
        for u in partial_words(e):
            for v in partial_words(e):
                entry = literal_components(u,v)
                if entry is None:
                    continue
                p, exceptional = entry
                if 2*exceptional != len(p):
                    continue
                a, b = len(u[1]), len(v[1])
                j = sum(y>=a for y in p[:a])
                d = (e-a-b+j)//2
                if d <= D:
                    result[d,e,a,b] += 1
    return result


def bounds_at(t, R, D, gf, shapes, allowed):
    """Bound each height pair before adding the pairs.

    Removing all stretched components leaves a minimum-length half. Its letters,
    the exceptional half and their interleaving determine the original word.
    The two side counts at each height pair may be multiplied for an upper
    bound, since this drops the requirement that the join be connected.
    """
    left = [Counter() for _ in range(D+1)]
    for (d, e, h, k), ways in shapes.items():
        if e>t:
            continue
        for (a,b,j), ordinary in gf[t-e].items():
            left[d][a+h,b+k] += comb(t,e)*ways*ordinary
    result=[]
    cumulative=Counter()
    for terms in left:
        cumulative.update(terms)
        result.append(sum(min(value,allowed[t][a,b])*allowed[R][a,b] for (a,b),value in cumulative.items()))
    return result


def run():
    sizes=((30,30),(30,37),(32,32),(32,40),(39,39),(39,47))
    limit=max(t for _,t in sizes)
    gf=gf_layers(limit,False)
    shapes=exceptional_shapes(2)
    baseline=single_factor_opposite(limit)
    patterns,distribution=forbidden_words(8)
    machine=automaton(patterns)
    allowed=possible_opposite(limit,machine)
    assert all(all(v<=baseline[t][a,b] for (a,b),v in counter.items()) for t,counter in enumerate(allowed))
    assert {x:v for x,v in shapes.items() if x[0]==1}=={(1,3,1,1):2,(1,4,2,0):3,(1,4,0,2):3}
    small=[]
    old_small=[]
    for t in range(1,9):
        hist=Counter();independent_allowed=Counter()
        for u in partial_words(t):
            for v in partial_words(t):
                # For at most eight letters there is a forbidden factor
                # precisely when the complete interval contains a closed cycle.
                upper,lower=u[0],v[0]
                w=tuple(2*int(0<=upper[i]<i)+int(0<=lower[i]<i) for i in range(t))
                if not internal_cycle(w):independent_allowed[len(u[1]),len(v[1])]+=1
                entry=literal_components(u,v)
                if entry is None:
                    continue
                p,_=entry
                a,b=len(u[1]),len(v[1])
                d=(t-a-b+sum(y>=a for y in p[:a]))//2
                hist[a,b,d]+=1
        assert independent_allowed==allowed[t]
        for R in range(t%2 or 2,9,2):
            predicted=bounds_at(t,R,2,gf,shapes,allowed)
            literal=[sum(v*allowed[R][a,b] for (a,b,d),v in hist.items() if d<=D) for D in range(3)]
            assert literal[0]==predicted[0] and all(x<=y for x,y in zip(literal,predicted))
            small.append([t,R,*literal])
            old_predicted=bounds_at(t,R,2,gf,shapes,baseline)
            old_literal=[sum(v*baseline[R][a,b] for (a,b,d),v in hist.items() if d<=D) for D in range(3)]
            assert old_literal[0]==old_predicted[0] and all(x<=y for x,y in zip(old_literal,old_predicted))
            old_small.append([t,R,*old_literal])
    # Exact terms of A005315; M_39 >= M_38 M_1 by the concatenation of closed meanders.
    denominators={30:1279842888607976889797535044,32:154967080111946539652848747626,
                  39:295598081822804246724487760690994952}
    rows=[]
    for n,t in sizes:
        old=bounds_at(t,2*n-t,2,gf,shapes,baseline)
        new=bounds_at(t,2*n-t,2,gf,shapes,allowed)
        assert all(x<=y for x,y in zip(new,old))
        rows.append(dict(n=n,side_length=t,opposite=2*n-t,max_excess=[0,1,2],
                         numerator_upper=new,single_forbidden_factor_numerator_upper=old,
                         denominator_lower=denominators[n],denominator_exact=n!=39))
    return dict(convention="One specified side; ordinary unweighted connected closed meanders; cumulative d.",
                exceptional_shapes=[[d,e,h,k,v] for (d,e,h,k),v in sorted(shapes.items())],
                small_checks=small,single_forbidden_factor_small_checks=old_small,orders=rows,
                avoidance=dict(maximum_factor_length=8,minimal_factor_counts=distribution,
                               states=len(machine),both_heights_nonnegative_at_every_step=True),
                sources={p.name:sha256(p.read_bytes()).hexdigest() for p in
                         (Path(__file__).resolve(),Path(__file__).with_name("minimum_caps.py"),
                          Path(__file__).with_name("check_excess_caps.py"),Path(__file__).with_name("height_envelope.py"))})


if __name__=="__main__":
    opt=ArgumentParser()
    opt.add_argument("--write",type=Path)
    args=opt.parse_args()
    data=run()
    if args.write:
        pending=args.write.with_suffix(".pending")
        pending.write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")
        pending.replace(args.write)
    else:
        home=Path(__file__).resolve()
        file=home.parent/"certificates/cut-excess-bound.json"
        if not file.is_file():
            file=home.parents[3]/"output/research73/star/cut-excess-bound.json"
        assert json.loads(file.read_bytes())==data
    for row in data["orders"]:
        print(row["n"],row["side_length"],[v/row["denominator_lower"] for v in row["numerator_upper"]])

