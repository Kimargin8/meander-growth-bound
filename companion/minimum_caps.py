"""Independent counts of shortest physical realisations of a pairing of the exposed ends."""
from argparse import ArgumentParser
from collections import Counter
from functools import cache
from itertools import product
from hashlib import sha256
from math import comb,factorial,prod
from pathlib import Path
import json

@cache
def matchings(k):
    if not k:return ((),)
    return tuple((2*i+1,*(x+1 for x in left),0,*(x+2*i+2 for x in right))
                 for i in range(k) for left in matchings(i) for right in matchings(k-i-1))

def split_pattern(order,p):
    indices={x:i for i,x in enumerate(order)}
    height=primitives=0
    for i,x in enumerate(order):
        height+=1 if indices[p[x]]>i else -1
        assert height>=0
        primitives+=height==0
    assert height==0
    return primitives

def prediction(a,p):
    # Labels go down the cut: outermost upper to innermost, then innermost lower to outermost.
    n=len(p); upper=list(range(a-1,-1,-1)); lower=list(range(a,n))
    first=[i for i,s in enumerate(upper) if p[s]>=a]; second=[i for i,s in enumerate(lower) if p[s]<a]
    assert [p[upper[i]] for i in first]==[lower[j] for j in second]
    cur_top=cur_low=0;orders=[]
    for i,j in zip([*first,len(upper)],[*second,len(lower)]):
        u,v=split_pattern(upper[cur_top:i],p),split_pattern(lower[cur_low:j],p)
        orders.append(comb(u+v,u))
        cur_top=i+1;cur_low=j+1
    return n-len(first),prod(orders)

@cache
def partial_words(t):
    out=[]
    for letters in product((False,True),repeat=t):
        arcs=[-1]*t; stack=[]
        for i,up in enumerate(letters):
            if up:stack.append(i)
            elif stack:
                j=stack.pop();arcs[j]=i;arcs[i]=j
            else:break
        else:out.append((tuple(arcs),tuple(stack)))
    return out

def trace_pairs(up,low):
    top,T=up;bottom,B=low
    stubs=[*( (s,0) for s in T),*( (s,1) for s in reversed(B))]
    if not stubs:return None
    where={s:i for i,s in enumerate(stubs)}
    mates=[-1]*len(stubs)
    seen=set()
    for i,(s,color) in enumerate(stubs):
        if mates[i]!=-1:continue
        while True:
            assert s not in seen
            seen.add(s); color=1-color
            nxt=(top,bottom)[color][s]
            if nxt<0:
                j=where[s,color];mates[i]=j;mates[j]=i;break
            s=nxt
    return (len(T),tuple(mates)) if len(seen)==len(top) else None

def one_cycle(p,q):
    if not p:return False
    s=0
    for k in range(len(p)):
        s=(p,q)[k%2][s]
        if s==0:return k+1==len(p)
    return False

def catalan(k):
    return comb(2*k,k)//(k+1)


def catalan_power(k,p):
    """Coefficient of x^k in C(x)^p."""
    if k<0:return 0
    if not p:return int(k==0)
    return p*comb(2*k+p,k)//(2*k+p)


def primitive_power(k,p):
    """p primitive Dyck factors of total semilength k."""
    return catalan_power(k-p,p)


def gf_layers(limit,canonical):
    # H is an arbitrary word of through, upper primitive and lower primitive
    # atoms. For S, each gap between through atoms is uniquely an upper
    # Dyck word followed by a lower Dyck word. The flag encodes that order.
    work=[Counter() for _ in range(limit+1)]
    work[0][0,0,0,False]=1
    for t in range(limit):
        for (a,b,j,lower_started),v in work[t].items():
            work[t+1][a+1,b+1,j+1,False]+=v
            for k in range(1,(limit-t)//2+1):
                weight=v*catalan(k-1)
                if not lower_started:
                    work[t+2*k][a+2*k,b,j,False]+=weight
                work[t+2*k][a,b+2*k,j,canonical]+=weight
    layers=[]
    for layer in work:
        out=Counter()
        for (a,b,j,_),v in layer.items():out[a,b,j]+=v
        layers.append(out)
    return layers


def closed_coefficients(t,a,b,j):
    # [w^j]H=(zuv)^j/(1-P(zu)-P(zv))^(j+1).
    # [w^j]S=(zuv)^j (C((zu)^2) C((zv)^2))^(j+1).
    if j>min(a,b) or t!=a+b-j or (a-j)%2 or (b-j)%2:return 0,0
    k=(a-j)//2;l=(b-j)//2
    types=catalan_power(k,j+1)*catalan_power(l,j+1)
    schedules=0
    for p in range(k+1):
        for q in range(l+1):
            orders=comb(p+q+j,j)*comb(p+q,p)
            schedules+=orders*primitive_power(k,p)*primitive_power(l,q)
    return schedules,types


def run():
    limit=11
    schedule_layers=gf_layers(limit,False)
    type_layers=gf_layers(limit,True)
    example=None
    observations=[]
    graded=[]
    total_checks=Counter()
    previous_even=1
    for t in range(limit+1):
        if t==0:
            hist=Counter()
        else:
            words=partial_words(t)
            hist=Counter()
            for u in words:
                for v in words:
                    tag=trace_pairs(u,v)
                    if tag:hist[tag]+=1
        total_checks["nonempty_histogram_entries"]+=len(hist)
        minima=legal=maxmult=0
        observed_schedules=Counter()
        observed_types=Counter()
        if t==0:
            # Both generating functions include their empty constant term.
            observed_schedules[0,0,0]=observed_types[0,0,0]=1
        seen_positive=0
        shortest_above_twelve=0
        max_all=0
        for a in range(t%2,t+1,2):
            for b in range(t%2,t+1,2):
                if not a+b:continue
                for m in matchings((a+b)//2):
                    L,count=prediction(a,m)
                    found=hist[a,m]
                    total_checks["reachability_equalities"]+=1
                    assert bool(found)==(t>=L),(a,b,t,L,m)
                    seen_positive+=bool(found)
                    if a+b<=12:legal+=bool(found)
                    if t!=L:continue
                    assert found==count,(a,b,t,count,found,m)
                    total_checks["shortest_multiplicity_equalities"]+=1
                    if a+b>12:shortest_above_twelve+=1
                    else:
                        maxmult=max(maxmult,count)
                        minima+=1
                    max_all=max(max_all,count)
                    j=sum(partner>=a for partner in m[:a])
                    assert t==a+b-j
                    observed_schedules[a,b,j]+=found
                    observed_types[a,b,j]+=1
                    if example is None and count>=6 and t==L and a+b<=12:
                        example=dict(a=a,b=b,t=t,pairing=m,minimum_schedules=count)
        # No observed nonplanar or otherwise unaccounted signature was skipped.
        assert seen_positive==len(hist)
        if t:
            observations.append(dict(t=t,word_pairs=len(words)**2,positive_types_up_to_twelve_ends=legal,
                                     shortest_types=minima,largest_shortest_multiplicity=maxmult))
        H=observed_schedules
        S=observed_types
        assert H==schedule_layers[t] and S==type_layers[t]
        positive_grades=[]
        checked=0
        for a in range(t+1):
            for b in range(t+1):
                for j in range(min(a,b)+1):
                    key=a,b,j
                    h,s=closed_coefficients(t,a,b,j)
                    assert H[key]==schedule_layers[t][key]==h
                    assert S[key]==type_layers[t][key]==s
                    if h or s:positive_grades.append([a,b,j,h,s])
                    checked+=1
        total_checks["four_degree_tuples_including_zero_coefficients"]+=checked
        total_checks["four_degree_nonzero_tuples"]+=len(positive_grades)
        total_h=sum(H.values())
        by_through=[]
        for j in range(t+1):
            h=sum(v for (a,b,number),v in H.items() if number==j)
            s=sum(v for (a,b,number),v in S.items() if number==j)
            if (t-j)%2:
                assert h==s==0
            else:
                k=(t-j)//2
                # [w^j]H(z,1,1,w)=z^j(1-4z²)^(-(j+1)/2).
                numerator=2**k*prod(j+1+2*i for i in range(k))
                coefficient,remainder=divmod(numerator,factorial(k))
                assert remainder==0 and h==coefficient
                assert s==catalan_power(k,2*(j+1))
                if t%2:
                    m=t//2;d=(j-1)//2
                    assert h==comb(m,d)*4**(m-d)
                    total_checks["odd_binomial_coefficients"]+=1
                by_through.append([j,h,s])
            total_checks["through_coefficients_per_series"]+=1
        if t%2:
            assert total_h==5**(t//2)
        elif t:
            previous_even=5*previous_even-2*catalan(t//2-1)
            assert total_h==previous_even
        else:assert total_h==1
        graded.append(dict(t=t,shortest_realisations=total_h,shortest_types=sum(S.values()),
                           shortest_types_above_twelve_ends=shortest_above_twelve,
                           largest_shortest_multiplicity=max_all if t else 1,
                           reachable_nonempty_types=len(hist),tested_degree_tuples=checked,
                           grades=positive_grades,through_profiles=by_through))
    # A first nonzero interaction can require the longest allowed physical test.
    a=b=2
    states=[(m,prediction(a,m)[0]) for m in matchings(2)]
    interior=[m for m,L in states if L<4]; edge=[m for m,L in states if L==4]
    assert len(interior)==len(edge)==1
    assert not one_cycle(interior[0],interior[0]) and one_cycle(interior[0],edge[0])
    ranks=dict(a=a,b=b,side_length=4,interior_pairing=interior[0],longest_pairing=edge[0],
               strict_interior_matrix=[[0]],interior_to_longest=[[1]])
    return dict(observations=observations,nontrivial_interleaving=example,strict_cut_example=ranks,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                all_grades=dict(columns=["upper_height","lower_height","through_pairs","shortest_realisations",
                                        "shortest_types"],through_columns=["through_pairs","shortest_realisations",
                                                                          "shortest_types"],
                                checks=dict(total_checks),layers=graded))

if __name__=="__main__":
    arg=ArgumentParser();arg.add_argument("--write",type=Path);args=arg.parse_args()
    data=run()
    if args.write:args.write.write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")
    print(sum(r["word_pairs"] for r in data["observations"]),data["nontrivial_interleaving"],data["all_grades"]["checks"])
