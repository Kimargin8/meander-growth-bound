"""An exact automaton certificate for the maximum cut height of closed meanders.

The ensemble is the ordinary, unit-weight enumeration by pairs of Dyck
words with fixed bank labels. Letters 0,1,2,3 represent UU, UD, DU, DD
in a scan from an outside end towards an interior cut.

A closed component inside a proper interval would be a second component
of the full meander. Local matched pairs can be computed with empty
stacks and unmatched pops ignored: a pair found this way persists for
every choice of earlier symbols. Thus both half-words at an interior cut
avoid every short factor containing an internal cycle. Reflecting the right
half reverses the order and interchanges U and D on both banks.

Write H=a+b for the height sum at the cut and Q>1 for the marking. The
letter weights Q,1,1,1/Q give each half the product Q^(H/2). We drop the
Dyck condition, connectivity and matching across the cut when upper-bounding
the two sums with the avoidance automaton. For its nonnegative weighted
matrix A, a positive vector v with Av <= r*v bounds the length-s weighted
sum, started at the empty state, by (v[0]/min(v))*r^s. Consequently the
number of meanders with H >= c*n at a given cut is at most
(v[0]/min(v))^2 * (r^2/Q^c)^n. There are 2*n-1 interior cuts; endpoints
have height zero. The n-th root of the probability of the union has upper
limit at most r^2/(mu*Q^c), using lim M_n^(1/n) >= mu = 11.380 as proved
by Albert and Paterson, JCTA 112 (2005), 250–262. Their unmarked language
method is a predecessor of this height-marked estimate.

Outside that exceptional set, at every cut any half of length ell >= n
has excess d=(ell-(a+b-j))/2, with j through pairs. Since j>=0 and
a+b<c*n<=c*ell, it follows that d/ell>(1-c)/2. At equal lengths both
sides obey the conclusion. These are counting estimates, with all
thresholds compared using integer powers rather than rounded logarithms.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import ceil
from pathlib import Path
import json


def internal_cycle(word):
    parents=list(range(len(word)))
    def root(i):
        while parents[i]!=i:
            parents[i]=parents[parents[i]]
            i=parents[i]
        return i
    for mask in (2,1):
        stack=[]
        for i,x in enumerate(word):
            if not x&mask:stack.append(i)
            elif stack:
                j=stack.pop()
                u,v=root(i),root(j)
                if u==v:return True
                parents[u]=v
    return False


def forbidden_words(limit):
    # Any proper factor lies entirely in the word without its first or
    # last letter. An already closed component persists on extension.
    previous=set();patterns=[];counts=[]
    for n in range(1,limit+1):
        now=set();size=0
        for letters in product(range(4),repeat=n):
            word=bytes(letters)
            inherited=word[1:] in previous or word[:-1] in previous
            if inherited or internal_cycle(word):
                now.add(word)
                if not inherited:
                    patterns.append(word);size+=1
        counts.append(size);previous=now
    known=set(patterns)
    assert all(bytes(3-x for x in reversed(w)) in known for w in patterns)
    return patterns,counts


def automaton(patterns):
    bad=set(patterns)
    prefixes=sorted({w[:i] for w in patterns for i in range(len(w))},key=lambda x:(len(x),x))
    ids={w:i for i,w in enumerate(prefixes)}
    transitions=[]
    for w in prefixes:
        row=[]
        for letter in range(4):
            s=w+bytes([letter])
            if any(s[i:] in bad for i in range(len(s))):
                row.append(-1)
                continue
            while s not in ids:s=s[1:]
            row.append(ids[s])
        transitions.append(row)
    return transitions


def support(rows,Q,powers):
    # Q=A/B. The integer matrix is AB times the weighted adjacency.
    a,b=Q.numerator,Q.denominator
    letters=[a*a,a*b,a*b,b*b]
    edges=[[(letter_weight,j) for letter_weight,j in zip(letters,r) if j>=0] for r in rows]
    def step(vector):
        return [sum(c*vector[j] for c,j in row) for row in edges]
    v=[1]*len(rows)
    for _ in range(powers):v=step(v)
    following=step(v)
    assert all(v)
    ratio=max(Fraction(x,y*a*b) for x,y in zip(following,v))
    # Round the super-eigenvalue upwards so the public comparison uses
    # short integers. The original positive vector is reconstructible.
    scale=10**9
    upper=Fraction(ceil(ratio*scale),scale)
    assert all(x*upper.denominator<=y*upper.numerator*a*b for x,y in zip(following,v))
    return upper,ratio,Fraction(v[0],min(v))


def critical_exponent(rho,Q,mu):
    # Seek the smallest per-mille height threshold with strict exponential decay.
    return next(Fraction(i,1000) for i in range(1,2001)
                if (rho*rho/mu)**1000<Q**i)


def rate_ceiling(rho,Q,c,mu):
    # Strict upward rounding to the next millionth, using only integer powers.
    exponent=c.denominator
    bound=(rho*rho/mu)**exponent/Q**c.numerator
    lo,hi=0,2_000_000
    while hi-lo>1:
        mid=(hi+lo)//2
        if Fraction(mid,1_000_000)**exponent>bound:hi=mid
        else:lo=mid
    result=Fraction(hi,1_000_000)
    assert result**exponent>bound
    return result


def prefix_scan(rows,word):
    s=0
    for x in word:
        s=rows[s][x]
        if s<0:return False
    return True


def dyck(n):
    def rec(p,h):
        if len(p)==2*n:
            if h==0:yield p
            return
        if h<2*n-len(p):yield from rec(p+(0,),h+1)
        if h:yield from rec(p+(1,),h-1)
    return list(rec((),0))


def ordinary_checks(rows,limit):
    totals=[]
    for n in range(1,limit+1):
        words=dyck(n)
        def pairing(word):
            q=[-1]*len(word);stack=[]
            for i,x in enumerate(word):
                if not x:stack.append(i)
                else:j=stack.pop();q[j]=i;q[i]=j
            return q
        systems=[(w,pairing(w)) for w in words]
        connected=0;maxima=Counter()
        for top,u in systems:
            for low,v in systems:
                i=0
                for k in range(2,2*n+1,2):
                    i=v[u[i]]
                    if not i:break
                if k!=2*n:continue
                connected+=1
                w=bytes(2*x+y for x,y in zip(top,low))
                # Both longest proper half-words must survive. Every shorter
                # cut is a prefix of one of these two scans.
                assert prefix_scan(rows,w[:-1]) and prefix_scan(rows,bytes(3-x for x in w[:0:-1]))
                h=hi=0
                for e in w:
                    h+=2*(e==0)-2*(e==3);hi=max(h,hi)
                assert h==0
                maxima[hi]+=1
        totals.append(dict(n=n,M=connected,max_height=dict(sorted(maxima.items()))))
    assert [r["M"] for r in totals]==[1,2,8,42,262,1828][:limit]
    return totals


def parse_options():
    p=ArgumentParser(description=__doc__)
    p.add_argument("--write",type=Path)
    return p.parse_args()


if __name__=="__main__":
    opts=parse_options()
    patterns,distribution=forbidden_words(8)
    mu=Fraction(569,50)
    results=[]
    for k,Q in ((2,Fraction(37,20)),(4,Fraction(7,4)),(6,Fraction(17,10)),(8,Fraction(17,10))):
        subset=[w for w in patterns if len(w)<=k]
        rows=automaton(subset)
        rho,raw,constant=support(rows,Q,48)
        c=critical_exponent(rho,Q,mu)
        bounds=[]
        for threshold in sorted({c,Fraction(63,100),Fraction(16,25),Fraction(17,25)}):
            if (rho*rho/mu)**threshold.denominator<Q**threshold.numerator:
                bounds.append(dict(height=str(threshold),longer_half_excess=str((1-threshold)/2),
                                   exceptional_fraction_nth_root_upper=str(rate_ceiling(rho,Q,threshold,mu))))
        results.append(dict(maximum_word_length=k,patterns=len(subset),automaton_states=len(rows),
                            Q=str(Q),iterations=48,rho_upper=str(rho),sharp_ratio_before_rounding=str(raw),
                            one_half_prefactor=str(constant),height_threshold=str(c),
                            longer_half_excess=str((1-c)/2),tradeoffs=bounds))
        print(k,len(subset),len(rows),float(Q),float(rho),float(c),float((1-c)/2),flush=True)
    # Compare the shortest-suffix transition automaton to a direct graph test,
    # including words whose central cycle matches neither bank at its ends.
    words_checked=0
    for size in range(1,9):
        for w in product(range(4),repeat=size):
            assert prefix_scan(rows,w)==(not internal_cycle(w))
            words_checked+=1
    output=dict(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                alphabet=["UU","UD","DU","DD"],mu_lower="569/50",
                mu_source="Albert and Paterson, Bounds for the growth rate of meander numbers, Theorem 1.1",
                mu_doi="10.1016/j.jcta.2005.02.006",factor_counts_by_length=distribution,
                minimal_forbidden=["".join(map(str,w)) for w in patterns],certificates=results,
                complete_word_checks=words_checked,small_closed=ordinary_checks(rows,6))
    if opts.write:
        tmp=opts.write.with_suffix(".pending")
        tmp.write_text(json.dumps(output,indent=2)+"\n",encoding="utf8");tmp.replace(opts.write)
    else:
        script=Path(__file__).resolve()
        path=script.parent/"certificates/height-envelope.json"
        if not path.is_file():path=script.parents[3]/"output/research73/star/height-envelope.json"
        assert json.loads(json.dumps(output))==json.loads(path.read_bytes())
