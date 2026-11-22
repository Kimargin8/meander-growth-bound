"""Certify two central-cut tail bounds for uniform connected, labelled closed meanders.

Write a,b for the two common endpoint heights. For the left or right half let
j count opposite-bank paths and put d=(n-a-b+j)/2. A completion of a fixed
half needs two length-n walks from empty stacks to these same heights when
read from the opposite outer end. At tilt s in (0,1), each bank has at most
(s+1/s)^n s^a (respectively with b) words: take its coefficient in the
binomial expansion. Ignoring positivity and connectedness increases the bound.

The generating function of shortest halves weighted by s per exposed end is
1/(sqrt(1-4s^2 z^2)-s^2 z), with coefficient growth bounded by C_s h^n,
h=s sqrt(4+s^2). This follows by rationalizing the denominator; the square-root
singularities lie outside the circle of the two simple poles. Removing all
oversized path components removes e<=4d positions, including two external
endpoints per removed component. Its actual weight s^(a+b) is therefore at
most the weight of its surviving shortest half. Positions and their two
binary steps recover the original words. For d<=n/100 on one half this gives
the upper n-th root B=(s+1/s)^2 * [4 h^24 * 25^25/24^24]^(1/25).
The exponent in the binomial tail increases up to p=1/25 for the stated tilt.
The union over left and right costs a factor of at most two, independent of n.

For a sharper deletion bound, let D=floor(n/100) and K be the net height
in the e erased upper-lower steps. Each of the m removed open paths has two
exposed ends, so K=2m. Each path has minimum length at most two, and all
the surplus is at most 2D. Hence K>=e-2D. For y>s the sum of the weights
s^K over allowed erased binary patterns is at most
(s/y)^(e-2D) (y+1/y)^(2e). After choosing the e positions, the remaining
length n-e contributes at most C h^(n-e). For s^2=5/16, y^2=5/3, the
binomial-entropy bound is increasing for 0<=e/n<=1/25, as checked below.
Its maximum at e/n=1/25 gives the sharper absolute base 7.841828909.

If min(d_left,d_right)<=0.09n, the shared heights satisfy a+b>=0.82n.
In each of the two half-word pairs the number of down steps is n-(a+b)/2,
so it is at most (59/200)(2n). Counting arbitrary binary words in both halves
gives growth at most exp(4 H(59/200)), with H the binary entropy using natural
logs. Dividing the two absolute growth bounds by mu>=11.380 gives the claimed
probability bounds. As a second check, the seam splice gives M_(k+m)>=M_k M_m:
exchange the upper arcs touching the last point of the left and the first
point of the right meander, joining those two points and their former partners.
The unique outer crossing arc and the new seam arc recover both inputs.
For any of the tail bounds, a simultaneous cyclic rotation of the labelled
sites in both matchings is a bijection. Allowing all possible origins of
half-length windows increases the bound by a factor of at most 2n.

Both bounds concern the limsup of n-th roots; a displayed root bound is not a
finite-n probability estimate. The certificate uses exact rational powers.
The independent small calculation tries pairs of full noncrossing matchings,
tests a single alternating cycle and traces both half geometries.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DEST=ROOT/"output/research73/star/central-excess.json"
PAPER="https://cs.otago.ac.nz/research/publications/oucs-2004-08.pdf"
JOURNAL="https://doi.org/10.1016/j.jcta.2005.02.006"


def str_frac(f):
    return [str(f.numerator),str(f.denominator)]


def bracket(powered,exponent,lower,upper):
    assert Fraction(lower)**exponent<powered<Fraction(upper)**exponent
    return dict(exponent=exponent,lower=lower,upper=upper)


def exact_bounds(m32):
    x=Fraction(1,3)  # s^2
    ap=Fraction(11380,1000)
    p=Fraction(1,25)
    h2=x*(4+x)
    pref=(1+x)**2/x
    assert h2==Fraction(13,9) and pref==Fraction(16,3)
    # d/dp [H(p)+p log(4/h)] is positive if 4(1-p)>hp.
    assert (4*(1-p))**2>h2*p*p
    b25=pref**25*4*h2**12*25**25/24**24
    entropy50=Fraction(200**200,59**59*141**141)
    entropy500=Fraction(2000**2000,593**593*1407**1407)
    sharp_x,y2=Fraction(5,16),Fraction(5,3)
    sharp_h2=sharp_x*(4+sharp_x)
    sharp_pref=(1+sharp_x)**2/sharp_x
    letter2=(sharp_x/y2)*((1+y2)**2/y2)**2
    assert letter2*(1-p)**2>sharp_h2*p*p
    sharp100=sharp_pref**100*sharp_h2**48*Fraction(25**25,24**24)**4*(y2/sharp_x)*letter2**2
    return dict(tilt_squared=str_frac(x),weighted_shortest_base_squared=str_frac(h2),
                completion_base=str_frac(pref),absolute_1_percent_25th_power=str_frac(b25),
                absolute_9_percent_50th_power=str_frac(entropy50),
                sharper_tilts_squared=[str_frac(sharp_x),str_frac(y2)],
                sharper_1_percent_100th_power=str_frac(sharp100),
                bounds=dict(absolute_1_percent=bracket(b25,25,"7.955633151","7.955633152"),
                            one_percent_AP=bracket(b25/ap**25,25,"0.699089","0.699090"),
                            one_percent_M32=bracket(b25**32/Fraction(m32)**25,800,"0.973822","0.973823"),
                            sharper_absolute_1_percent=bracket(sharp100,100,"7.841828908","7.841828909"),
                            sharper_one_percent_AP=bracket(sharp100/ap**100,100,"0.689088","0.689089"),
                            sharper_one_percent_M32=bracket(sharp100**8/Fraction(m32)**25,800,"0.959892","0.959893"),
                            absolute_9_percent=bracket(entropy50,50,"11.316619","11.316620"),
                            nine_percent_AP=bracket(entropy50/ap**50,50,"0.994430","0.994431"),
                            nine_point_three_percent_AP=bracket(entropy500/ap**500,500,"0.999620","0.999621")),
                denominator=dict(AP="11.380",normalization="lim M_n^(1/n), 2n labelled sites",
                                 theorem="Theorem 1.1",author_copy=PAPER,journal=JOURNAL,
                                 m32=str(m32)))


@lru_cache(None)
def matchings(n):
    if not n:
        return ((),)
    result=[]
    for k in range(n):
        for inside in matchings(k):
            first=(2*k+1,)+tuple(v+1 for v in inside)+(0,)
            for tail in matchings(n-k-1):
                result.append(first+tuple(v+2*k+2 for v in tail))
    assert len(result)==comb(2*n,n)//(n+1)
    return tuple(result)


def signed_words(m):
    n=len(m)//2
    return (sum((m[i]>i)<<i for i in range(n)),
            sum((m[2*n-1-i]<2*n-1-i)<<i for i in range(n)))


def prefix_graph(n,upper,lower):
    links=[[None]*n for _ in range(2)]
    ends=[]
    for side,word in enumerate((upper,lower)):
        stack=[]
        for i in range(n):
            if word&(1<<i):
                stack.append(i)
            else:
                back=stack.pop()
                links[side][i],links[side][back]=back,i
        ends.append(tuple(stack))
    a,b=map(len,ends)
    visited=set()
    through=0
    for side in (0,1):
        for start in ends[side]:
            if start in visited:
                continue
            i=start
            visited.add(i)
            next_side=1-side
            while links[next_side][i] is not None:
                i=links[next_side][i]
                assert i not in visited
                visited.add(i)
                next_side=1-next_side
            through+=(next_side!=side)
    if len(visited)<n:
        return None
    excess=n-a-b+through
    down=2*n-upper.bit_count()-lower.bit_count()
    assert excess>=0 and excess%2==0 and 2*down==2*n-a-b
    assert a+b>=n-excess and 2*down<=n+excess
    return (a,b,through,excess//2,down)


def shortest_weight(n):
    values=[Fraction(1)]
    for t in range(1,n+1):
        value=values[t-1]/3
        for k in range(1,t//2+1):
            primitive=comb(2*k-2,k-1)//k
            value+=Fraction(2*primitive,3**k)*values[t-2*k]
        values.append(value)
    return values[n]


def check_order(n,known):
    mats=matchings(n)
    all_codes=[signed_words(m) for m in mats]
    codes=sorted({code for pair in all_codes for code in pair})
    physical={(x,y):prefix_graph(n,x,y) for x in codes for y in codes}
    shortest=Counter()
    for row in physical.values():
        if row is not None:
            a,b,j,d,_=row
            if d==0:
                shortest[(a+b)//2]+=1
    weight=sum(Fraction(v,3**e) for e,v in shortest.items())
    assert weight==shortest_weight(n)

    # Pair independent full diagrams. Exchange of the two banks preserves every
    # recorded statistic, including each side's excess. Off the diagonal the
    # two orderings are distinct. This avoids executing the same cycle twice.
    joint=Counter()
    min_hist=Counter()
    bright=Counter()
    for i,top in enumerate(mats):
        lx,rx=all_codes[i]
        for k in range(i+1):
            low=mats[k]
            u=low[top[0]]
            period=1
            while u:
                u=low[top[u]]
                period+=1
            if period!=n:
                continue
            ly,ry=all_codes[k]
            a,b,jl,dl,down_l=physical[lx,ly]
            aa,bb,jr,dr,down_r=physical[rx,ry]
            assert (a,b)==(aa,bb) and down_l==down_r
            least=min(dl,dr)
            assert a+b>=n-2*least and 2*down_l<=n+2*least
            if 100*least<=9*n:
                assert 50*(a+b)>=41*n and 100*down_l<=59*n
            multiple=1+(i!=k)
            joint[dl,dr]+=multiple
            min_hist[least]+=multiple
            # Each alpha appears in this count with its correct upper-lower ordering.
            bright[lx,ly]+=1
            if i!=k:
                bright[ly,lx]+=1
    assert sum(min_hist.values())==sum(joint.values())==sum(bright.values())==known
    assert all(joint[r,l]==v for (l,r),v in joint.items())
    for (x,y),number in bright.items():
        a,b,_,_,_=physical[x,y]
        ub=comb(n,(n-a)//2)*comb(n,(n-b)//2)
        assert number<=ub and ub*3**(n+(a+b)//2)<=16**n
    return dict(n=n,closed=known,ordered_diagram_pairs=len(mats)**2,
                min_excess={str(i):min_hist[i] for i in sorted(min_hist)},
                left_excess={str(i):sum(v for (x,y),v in joint.items() if x==i) for i in sorted({x for x,y in joint})},
                joint_excess=[[l,r,v] for (l,r),v in sorted(joint.items())],
                any_1_percent=sum(v for d,v in min_hist.items() if 100*d<=n),
                any_9_percent=sum(v for d,v in min_hist.items() if 100*d<=9*n),
                literal_weighted_shortest=str_frac(weight))


def check():
    if not __debug__:
        raise ValueError("Run with Python assertions enabled")
    ref=ROOT/"references/oeis/2026-11-15/b005315.txt"
    known={int(i):int(v) for i,v in (s.split() for s in ref.read_text().splitlines() if s.strip())}
    return dict(model="Uniform connected ordered pair of noncrossing perfect matchings on 2n labelled sites",
                event="min(d_left,d_right)<=delta*n at the central cut",
                asymptotic=exact_bounds(known[32]),literal=[check_order(n,known[n]) for n in range(1,10)],
                inputs={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),ref)})


if __name__=="__main__":
    result=check()
    if "--write" in sys.argv:
        DEST.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    else:
        assert json.loads(DEST.read_bytes())==result
    print([(r["n"],r["closed"],r["min_excess"]) for r in result["literal"]])
