"""Verify the exact certificate, the return bound, and the finite accuracy constants."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import re
from build import DEN, ROOT, TARGET, tables

OUT=ROOT/"output/research75/growth"
OLD=F(4581442652581,392062763189)
OLD_C=F(471986,10**12)
R=F(3898274927,10**9)
Q=F(159,100)


def endpoint_factor(m,i,j,x):
    """s^(|j-i|+1) times the corresponding entry of (sI-P)^(-1), x=s^2."""
    low,high=sorted((i,j))
    k=high-low+1
    det=[1,1]
    for p in range(2,m+1):
        det.append((x if p%2==0 else 1)*det[-1]-det[-2])
    exponent=k+(low%2)+(m-1-high)%2-m%2
    assert exponent%2==0
    return x**(exponent//2)*det[low]*det[m-1-high]/det[m]


def terminal_factor(upper_weight):
    x0=F(117,10)
    candidates=[1/(1+upper_weight)**6]
    # The parameter s^k G is decreasing in s^2 by its walk expansion.
    for a in range(1,33):
        for target in (a,a+1):
            lo,hi=max(0,target-4),min(34,target+4)
            for current in range(lo,hi+1):
                k=1+abs(current-target)
                scaled=endpoint_factor(hi-lo+1,current-lo,target-lo,x0)
                candidates.append(1/((1+upper_weight)**(6-k)*scaled))
    return min(candidates)


def finite_width(n,epsilon,constant,x):
    target=(R*R/x)**n * (2*(2*n-1)/(epsilon*constant))
    lo,hi=-1,n
    while hi-lo>1:
        k=(lo+hi)//2
        if Q**(2*k)>=target:
            hi=k
        else:
            lo=k
    return hi*2


def run():
    direct,walk,radius=tables()
    assert direct[16]==walk[16]==radius[16]
    stored=json.loads((OUT/"weights.json").read_bytes())
    assert stored["width32_all"]==direct[16] and stored["width32_radius4"]==radius[16]
    assert stored["denominator"]==DEN and F(stored["squared_parameter"])==TARGET
    one=(OUT/"iteration.txt").read_text()
    two=(OUT/"audit.txt").read_text()
    assert "width 32 cells 810375650 root 1 odd 386058577 even 424317073" in one
    assert "width 32 round 440 odd 386058577 even 424317073" in two
    assert [int(row) for row in re.findall(r"round (\d+) min ",one)]==list(range(360,441,5))
    x=F(re.search(r"round 440 min (\d+/\d+)",one)[1])
    assert x==F(re.search(r"\bmin (\d+/\d+)",two)[1]) and F(117,10)<x<TARGET
    r,V=map(int,re.search(r"root (\d+) largest (\d+)",two).groups())
    assert (R*R/x)**25<Q**14
    previous=ROOT/"output/research74/boundary"
    middle=F(re.search(r"\bmin (\d+/\d+)",(previous/"w32-audit.txt").read_text())[1])
    assert middle==F(5660040998411745,483083866680971)
    assert (R*R/middle)**25 > Q**14
    assert Q**3>4
    comp=(OUT/"compare.txt").read_text()
    oldroot,newroot=map(int,re.search(r"roots (\d+) (\d+)",comp).groups())
    g=F(re.search(r"rows 386058577 min (\d+/\d+)",comp)[1])
    old=(ROOT/"output/research73/growth/finite32.txt").read_text()
    assert f"round 300 odd 386058577 shift 4 root {oldroot} " in old
    assert newroot==r
    g40=F(re.search(r"domination 40 1 0 (\d+/\d+)",old)[1])
    exponent=int(re.search(r"level 40 exponent (\d+)",old)[1])
    beta=terminal_factor(F(max(direct[16]),DEN))
    gamma=2**exponent*g40*g
    finite=min(gamma*r*beta/x**41,OLD_C*(OLD/x)**43)
    lower=F(int(10**14*finite),10**14)
    assert lower<=finite and lower>0
    examples=[]
    for n in (39,100,1000,20000):
        eps=F(1,100)
        a=finite_width(n,eps,OLD_C,OLD)
        b=finite_width(n,eps,lower,x)
        examples.append(dict(n=n,epsilon=str(eps),published_width=a,new_width=b,combined=min(a,b)))
    checkpoints={}
    for path in (previous/"checkpoint-sha256.txt",OUT/"checkpoint.sha256"):
        for line in path.read_text().splitlines():
            digest,name=line.split()
            assert re.fullmatch(r"[0-9a-f]{64}",digest)
            checkpoints[Path(name).name]=digest
    assert set(checkpoints)=={"width32.local.checkpoint","width32-round360.bin","width32-round440.bin"}
    return dict(odd=386058577,even=424317073,round=440,x=str(x),root=r,maximum=V,
                checkpoints_sha256=checkpoints,
                parameter_upper=str(TARGET),edge_denominator=DEN,edge_numerators=direct[16],
                normalized_old_over_new=str(g),old_root_lower_multiplier=str(g40),
                old_denominator_exponent=exponent,auxiliary_radius=4,terminal_six_step_factor=str(beta),
                constant_for_all_orders=str(lower),rational_n_exponent="14/25",
                finite_width_examples=examples)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    receipt=OUT/"verified.json"
    if args.write:
        receipt.write_text(json.dumps(result,indent=2)+"\n")
    else:
        assert result==json.loads(receipt.read_bytes())
    print(json.dumps(result,indent=2))
