"""Audit the local integer transcripts or regenerate them using a C++17 compiler."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "output/research73/growth"
ORIGINAL = ROOT / "research/review57/strip/confined.cpp"
OEIS = ROOT / "references/oeis/2026-11-15/b005315.txt"
RECEIPT = OUT / "growth-bounds.json"
RUNS = {
    "w14": ("strip_bound.cpp", 14, 100),
    "w16": ("strip_bound.cpp", 16, 100),
    "w18": ("strip_bound.cpp", 18, 100),
    "w26": ("strip_bound.cpp", 26, 250),
    "f26": ("folded_bound.cpp", 26, 250),
    "f28": ("folded_bound.cpp", 28, 300),
    "f30": ("folded_bound.cpp", 30, 320),
    "i32": ("implicit_bound.cpp", 32, 300),
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse(text, width):
    lines = text.splitlines()
    start = lines[0].split()
    assert int(start[1]) == width
    counts = {int(row.split()[1]): int(row.split()[2]) for row in lines if row.startswith("count ")}
    brackets = {}
    for row in lines:
        if row.startswith("round "):
            _, i, _, low, _, _, _, high, _, _ = row.split()
            lower, upper = Fraction(low), Fraction(high)
            assert 0 < lower <= upper <= 16
            brackets[int(i)] = [str(lower), str(upper)]
    assert brackets
    if start[0] == "W":
        states = sum(comb(2*k, k)//(k+1)*(2*k+1) for k in range(1, width//2+1))
        assert [int(start[i]) for i in (3,5,7,9)] == [states]*4
        geometry = {"states": states, "root_communicating_class": states}
    else:
        assert start[0] == "width"
        cells = sum(comb(2*k, k) for k in range(1, width//2+1))
        odds = sum(comb(2*k, k)//(k+1)*((k+1)//2) for k in range(1, width//2+1))
        assert int(start[3]) == cells and int(start[5]) == 1
        assert int(start[7]) == odds and int(start[9]) == cells-odds
        returned = [int(row.split()[1]) for row in lines if row.startswith("all_odd_return_in")]
        assert len(returned) == 1 and returned[0] > 0 and not returned[0]%2
        geometry = {"states": cells, "odd_states": odds, "common_return_moves": returned[0]}
    return geometry, counts, brackets


@lru_cache(maxsize=None)
def arches(n):
    if not n:
        return ((),)
    result = []
    for k in range(n):
        for inner in arches(k):
            left = (2*k+1,) + tuple(x+1 for x in inner) + (0,)
            for other in arches(n-1-k):
                result.append(left + tuple(x+2*k+2 for x in other))
    return tuple(result)


def literal(n):
    mats = arches(n)
    heights = []
    for p in mats:
        h = 0
        row = []
        for i, j in enumerate(p):
            h += 1 if i < j else -1
            row.append(h)
        heights.append(row)
    hist = Counter()
    for i, p in enumerate(mats):
        for j in range(i+1):
            q = mats[j]
            v, visits = q[p[0]], 1
            while v:
                v, visits = q[p[v]], visits+1
            if visits == n:
                width = max(a+b for a,b in zip(heights[i],heights[j]))
                hist[width] += 1+(i!=j)
    return {w: sum(count for k,count in hist.items() if k<=w) for w in range(2, 2*n+1, 2)}


def compile_one(cpp, build, cxx):
    binary = build / cpp.removesuffix(".cpp")
    if not binary.exists():
        subprocess.run([cxx, "-std=c++17", "-O3", "-fopenmp", str(HERE/cpp),
                        "-o", str(binary)], check=True)
    return binary


def edge_summary(text):
    totals = Counter()
    for k,line in enumerate(text.splitlines()[:-1],1):
        fields=line.split()
        width,matchings,shifts,opens,closes,forbidden=map(int,fields[1::2])
        catalan=comb(2*k,k)//(k+1)
        previous=comb(2*k-2,k-1)//k
        assert width==2*k and matchings==catalan and shifts==(2*k+1)*catalan
        assert opens==(k+1)*catalan*(k<16)
        assert forbidden==k*previous and closes==k*(catalan-previous)
        totals.update(dict(states=catalan*(k+1),matchings=matchings,shifts=shifts,
                           opens=opens,closes=closes,early_closures=forbidden))
    last=text.splitlines()[-1].split()
    assert last[0]=="TOTAL" and len(text.splitlines())==17
    assert dict(totals)==dict(zip(last[1::2],map(int,last[2::2])))
    return dict(totals)


def finite_data(text):
    lines=text.splitlines()
    label=lines[0].split()
    info=dict(zip(label[1::2], map(int,label[2::2])))
    levels={}
    comparisons={}
    candidates={}
    for row in lines[1:]:
        v=row.split()
        if v[0]=="level":
            levels[int(v[1])]=[int(v[3]),int(v[5]),int(v[7])]
        elif v[0]=="domination":
            m,q,p=map(int,v[1:4])
            comparisons[q,p,m]=Fraction(v[4])
        elif v[0]=="candidate":
            q,p=map(int,v[1:3])
            candidates[q,p]=(Fraction(v[4]),Fraction(v[6]),int(v[8]),int(v[10]),int(v[12]))
        else:
            raise ValueError(row)
    return info,levels,comparisons,candidates


def floor_scientific(value):
    # Return a short, exact decimal lower bound to six significant digits.
    exponent=0
    while value<1:
        value*=10;exponent-=1
    while value>=10:
        value/=10;exponent+=1
    return {"integer":int(value*100000),"power_of_ten":exponent-5}


def finite_summary(text, lower, upper, known):
    info,levels,ratios,candidates=finite_data(text)
    assert info["width"]==32 and info["round"]==300 and info["odd"]==386058577
    assert 0<info["minimum"]<=info["root"]<=info["maximum"]<=2**52
    assert set(levels)==set(range(41)) and set(candidates)=={(1,0),(2,23),(5,58)}
    assert set(ratios)=={(q,p,m) for q,p in candidates for m in (15,40)}
    assert candidates[1,0]==(lower,upper,info["root"],info["minimum"],info["maximum"])
    for m,(exponent,vroot,largest) in levels.items():
        assert vroot<=largest and vroot>0 and exponent>=0
        if m<=15:
            assert exponent==0 and known[m+1]==vroot
        if m+1 in known:
            assert (vroot<<exponent)<=known[m+1]
    results={}
    for (q,p),(lo,hi,root,least,largest) in sorted(candidates.items()):
        assert 0<lo<=hi and 0<least<=root<=largest
        raw=Fraction(root,largest)
        exact15=root*ratios[q,p,15]
        assert raw<=exact15<=levels[15][1]
        finite_rates={}
        for key,L in (("exact_L",lo),("L_rounded_down_2_places",Fraction(int(100*lo),100))):
            tail=(2**levels[40][0])*root*ratios[q,p,40]/L**41
            possibilities={str(m+1):Fraction(r<<e)/L**(m+1) for m,(e,r,_) in levels.items() if m<40}
            possibilities["tail_n_at_least_41"]=tail
            minimizer=min(possibilities,key=possibilities.get)
            c=possibilities[minimizer]
            finite_rates[key]={"L":str(L),"c_for_every_n_at_least_1":str(c),"short_lower_for_c":floor_scientific(c),
                               "minimum_at":minimizer, "c_above_implies_Mn_at_least":"c * L^n"}
        results[f"{q},{p}"]={"positive_witness":"v" if p==0 else f"{q}Av-{p}v",
                              "root":root,"minimum":least,"maximum":largest,
                              "lower_ratio":str(lo),"upper_ratio":str(hi),"root_over_max":str(raw),
                              "K15":str(exact15),"K15_statement":"M_n >= K15 * L^(n-16) for n>=16",
                              "root_vector_comparisons":{str(m):str(ratios[q,p,m]) for m in (15,40)},
                              "finite_uniform":finite_rates}
    best_lo=max(candidates,key=lambda t:candidates[t][0])
    best_hi=min(candidates,key=lambda t:candidates[t][1])
    return {"witness":info,"candidates":results,
            "strongest_lower":{"q,p":list(best_lo),"exact":str(candidates[best_lo][0])},
            "strongest_strip_upper":{"q,p":list(best_hi),"exact":str(candidates[best_hi][1])}}


def cross_check(plain,finite):
    """A separate implementation numbers its own words and recomputes the matrix action."""
    lines=plain.splitlines()
    assert len(lines)==4
    info,_,_,candidates=finite_data(finite)
    fields=lines[0].split()
    assert fields[0::2]==["width","round","odd","even","root","largest"]
    assert list(map(int,fields[1::2]))==[32,300,info["odd"],424317073,info["root"],info["maximum"]]
    for row in lines[1:]:
        data=row.split()
        if data[0]=="shift":
            t=Fraction(data[1])
            q,p=t.denominator,t.numerator
            data=data[2:]
        else:
            q,p=1,0
        assert data[0]=="min" and data[2]=="max"
        assert (Fraction(data[1]),Fraction(data[3]))==candidates[q,p][:2]
    return {"independent_program":"verify_vector.cpp","agreement":"all six exact ratios, root coordinate and maximum"}


def native_tests(cxx, threads, selected, small, edges, finite):
    env = {**os.environ, "OMP_NUM_THREADS": str(threads)}
    with TemporaryDirectory() as temporary:
        build = Path(temporary)
        if small:
            table = {n:literal(n) for n in range(1,9)}
            for width in range(2,15,2):
                for cpp in ("strip_bound.cpp","folded_bound.cpp","implicit_bound.cpp"):
                    binary = compile_one(cpp, build, cxx)
                    log = subprocess.check_output([binary,str(width),"10"], env=env,text=True)
                    _,counts,_ = parse(log,width)
                    for n,row in table.items():
                        assert counts[n] == row[min(width,2*n)]
            print("All three transfers match direct labelled diagrams, orders 1 through 8, widths 2 through 14.")
        for key in selected:
            cpp,width,rounds = RUNS[key]
            binary = compile_one(cpp, build, cxx)
            fresh = subprocess.check_output([binary,str(width),str(rounds)], env=env,text=True)
            expected = parse((OUT/(key+".txt")).read_text(),width)
            assert parse(fresh,width) == expected, key
            print(key, "all exact ratios, width-restricted counts and return tests reproduced")
        if edges:
            binary=compile_one("check_strip_edges.cpp",build,cxx)
            fresh=subprocess.check_output([binary,"32"],env=env,text=True)
            assert edge_summary(fresh)==edge_summary((OUT/"edges32.txt").read_text())
            print("All width-32 edges decoded as literal matchings.")
        if finite:
            binary=compile_one("implicit_bound.cpp",build,cxx)
            checkpoint=build/"odd.checkpoint"
            fresh=subprocess.check_output([binary,"32","300",checkpoint],env=env,text=True)
            assert parse(fresh,32)==parse((OUT/"i32.txt").read_text(),32)
            reader=compile_one("finite_bound.cpp",build,cxx)
            log=subprocess.check_output([reader,checkpoint,"40"],env=env,text=True)
            assert finite_data(log)==finite_data((OUT/"finite32.txt").read_text())
            print("Width-32 witness and finite lower estimates independently regenerated.")


def check():
    known = dict(tuple(map(int,s.split())) for s in OEIS.read_text().splitlines() if s.strip())
    trials = {}
    inputs = [Path(__file__).resolve(), ORIGINAL, OEIS,HERE/"check_strip_edges.cpp",
              HERE/"finite_bound.cpp",OUT/"edges32.txt",OUT/"finite32.txt",
              HERE/"verify_vector.cpp",OUT/"vector32-final.txt"]
    for key,(source,width,nsteps) in RUNS.items():
        path = OUT / (key+".txt")
        geometry,counts,brackets = parse(path.read_text(),width)
        assert max(brackets) == nsteps and set(counts) == set(range(1,13))
        for n,c in counts.items():
            # The two Dyck heights are <=n, so W >= 2n suffices for all pairs.
            if 2*n<=width:
                assert known[n] == c
        lo,hi = (Fraction(t) for t in brackets[nsteps])
        trials[key] = {**geometry,"width":width,"iterations":nsteps,"counts_1_to_12":list(counts.values()),
                       "lower_ratio":str(lo),"upper_ratio":str(hi),
                       "lower_truncated_6_places":int(lo*10**6),
                       "upper_ceiling_6_places":-((-hi.numerator*10**6)//hi.denominator)}
        inputs.extend([path,HERE/source])
    assert trials["w26"]["counts_1_to_12"] == trials["f26"]["counts_1_to_12"]
    assert Fraction(trials["w26"]["upper_ratio"]) >= Fraction(trials["f26"]["lower_ratio"])
    assert Fraction(trials["f26"]["upper_ratio"]) >= Fraction(trials["w26"]["lower_ratio"])
    assert Fraction(trials["w26"]["lower_ratio"]) > Fraction(1144,100)
    assert Fraction(trials["f28"]["lower_ratio"]) > Fraction(1154,100)
    assert Fraction(trials["f30"]["lower_ratio"]) > Fraction(116,10)
    assert Fraction(trials["i32"]["lower_ratio"]) > Fraction(116,10)
    small = {str(n): {str(k):v for k,v in literal(n).items()} for n in range(1,9)}
    assert all(known[int(n)] == row[str(2*int(n))] for n,row in small.items())
    return {"model":"uniform labelled ordered connected meanders; at most W total upper and lower open ends",
            "normalization":"C_n(W)=(T^(2n-2))_oo; bounds are per n",
            "trials":trials, "direct_counts_upper_width_by_order":small,
            "all_edges_at_width_32":edge_summary((OUT/"edges32.txt").read_text()),
            "uniform_finite_bound":finite_summary((OUT/"finite32.txt").read_text(),
                                                   Fraction(trials["i32"]["lower_ratio"]),
                                                   Fraction(trials["i32"]["upper_ratio"]),known),
            "independent_vector_check":cross_check((OUT/"vector32-final.txt").read_text(),
                                                   (OUT/"finite32.txt").read_text()),
            "sources_and_transcripts":{p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(inputs))}}


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("enable Python assertions")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--small-native",action="store_true")
    ap.add_argument("--edges32",action="store_true")
    ap.add_argument("--finite32",action="store_true")
    ap.add_argument("--replay",action="append",choices=RUNS,default=[])
    ap.add_argument("--cxx",default="g++")
    ap.add_argument("--threads",type=int,default=4)
    opts = ap.parse_args()
    if opts.small_native or opts.replay or opts.edges32 or opts.finite32:
        native_tests(opts.cxx,opts.threads,opts.replay,opts.small_native,opts.edges32,opts.finite32)
    else:
        result=check()
        if opts.write:
            RECEIPT.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
        else:
            assert result == json.loads(RECEIPT.read_bytes())
        print("Exact transcripts and direct diagram checks verified.")
