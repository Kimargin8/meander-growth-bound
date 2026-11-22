"""Check rooted return witnesses with integers and physical small chord gluings."""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[2] / "output/research73/growth_return"
if (HERE / "certificates").is_dir():
    DATA = HERE / "certificates"


@lru_cache(None)
def systems(n):
    if n == 0:
        return ((),)
    result = []
    for i in range(n):
        partner = 2 * i + 1
        for left in systems(i):
            for right in systems(n - 1 - i):
                result.append((partner,) + tuple(1+j for j in left) + (0,)
                              + tuple(1+partner+j for j in right))
    return tuple(result)


def direct_histogram(n):
    """Follow alternating edges around point zero in a full planar diagram."""
    hist = Counter()
    choices = [(matching, tuple(1 if i < j else -1 for i, j in enumerate(matching)))
               for matching in systems(n)]
    for upper, u in choices:
        for lower, v in choices:
            point = 0
            for distance in range(1, 2*n+1):
                point = upper[point] if distance % 2 else lower[point]
                if point == 0:
                    break
            if distance == 2*n:
                h = peak = 0
                for a,b in zip(u,v):
                    h += a+b
                    peak = max(peak,h)
                hist[peak] += 1
    return hist


def parse(name):
    rows = [line.split() for line in (DATA / name).read_text().splitlines()]
    headers = [list(map(int, r[1:])) for r in rows if r[0] == "graph"]
    seeds = [int(r[1]) for r in rows if r[0] == "seed"]
    small = {int(r[1]):int(r[2]) for r in rows if r[0] == "small"}
    returns = [list(map(int,r[1:])) for r in rows if r[0] == "return"]
    assert len(headers) == len(seeds) == 1 and len(small) == 10 and returns
    assert all(a[0] < b[0] for a,b in zip(returns,returns[1:]))
    return headers[0],seeds[0],small,returns


def decimal_lower(m, x, exponent, seed):
    scale = 10**6
    lower, upper = 0, 16*scale
    # Returns the largest numerator at this decimal scale with strict inequality.
    target = (x << exponent) * scale**m
    while upper-lower > 1:
        mid = (lower+upper)//2
        if seed * mid**m < target:
            lower = mid
        else:
            upper = mid
    assert seed * lower**m < target <= seed * upper**m
    return lower


def certificate():
    base = {n:direct_histogram(n) for n in range(1,8)}
    counts = [sum(x.values()) for x in base.values()]
    assert counts == [1,2,8,42,262,1828,13820]
    hist = {}
    file_names = []
    def log(name):
        file_names.append(name)
        return parse(name)
    for width in (2,4,6,8,10,12):
        header,_,small,_ = log(f"small{width}.txt")
        assert header[0] == width
        expected = [sum(k for w,k in base[n].items() if w <= width) for n in base]
        assert expected == [small[n] for n in base]
        hist[width] = expected
    graph,seed,small,states = log("return26.txt")
    assert graph[0] == 26 and small == dict(enumerate(
        [1,2,8,42,262,1828,13820,110954,933458,8152860], start=1))
    scale = 10**6
    witnesses = []
    for m,x,e,maximum,cap_count in states:
        assert cap_count == 0 and 0 < x <= maximum <= 2**47
        if m in (1400,1800,1900,2500,3100):
            num = decimal_lower(m,x,e,seed)
            witnesses.append(dict(slots=m,root=x,binary_exponent=e,
                                  growth_strictly_greater_than=f"{num//scale}.{num%scale:06d}"))
    assert len(witnesses) == 5 and states[-1][0] == 3100
    fixed_graph,fixed_seed,fixed_small,returns = log("fixed26.txt")
    assert fixed_graph == graph and fixed_seed == seed and fixed_small == small
    m,x,e,maximum,cap_count = returns[-1]
    assert e == 0 and x > fixed_seed and maximum <= 2**50-1
    return dict(convention="n is half the number of crossings; a rooted return in m two-site slots has order m+1",
                width=26,even_states=graph[1],odd_states=graph[2],root_in_odd=graph[3],
                upper_and_lower_banks_are_ordered=True,reflection_quotient=False,scc_pruning=False,
                seed=seed,direct_width_checks=42,small_width_counts=hist,binary_witnesses=witnesses,
                fixed_witness=dict(numerator=57,denominator=5,slots=m,root=x,
                                   max_value=maximum,saturated_values=cap_count),
                sha256={**{name:sha256((DATA/name).read_bytes()).hexdigest() for name in file_names},
                        **{p:sha256((HERE/p).read_bytes()).hexdigest() for p in ("return.cpp","verify.py")}})


if __name__ == "__main__":
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--write",action="store_true")
    args = ap.parse_args()
    result = json.loads(json.dumps(certificate()))
    path = DATA / "verified.json"
    if args.write:
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
        tmp.replace(path)
    else:
        assert json.loads(path.read_bytes()) == result
    print(result["binary_witnesses"][-1],result["fixed_witness"])
