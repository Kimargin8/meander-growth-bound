"""Exact checks for a finite denominator and relative approximation at prescribed accuracy."""
from argparse import ArgumentParser
from fractions import Fraction
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/research74/boundary"
Q = Fraction(159,100)
R = Fraction(3898274927,10**9)
X = Fraction(5660040998411745,483083866680971)
OLD = Fraction(4581442652581,392062763189)
OLD_FACTOR = Fraction(471986,10**12)
sys.path.insert(0, str(ROOT/"research/review73/star"))


def width(n, eps, denominator):
    if n < 1 or not 0 < eps < 1:
        raise ValueError("Expected a positive order and an accuracy strictly between zero and one.")
    # The avoidance estimate uses only proper factors in both halves of the cut.
    # For small orders the full width gives the exact result automatically.
    if n < 16:
        return 2*n
    target = 2*(2*n-1)*R**(2*n)/(eps*denominator)
    left,right = -1,n
    while right-left>1:
        m=(left+right)//2
        if Q**(2*m) >= target:
            right=m
        else:
            left=m
    assert right==n or Q**(2*right)>=target
    assert right==0 or Q**(2*right-2)<target
    return 2*right


def check():
    from approximation import cycle_language
    from height_envelope import automaton, support
    from strip_error import connected_strip
    assert connected_strip(6,(1,2,3)) == {1:1,2:2,3:8}
    patterns, _ = cycle_language(10)
    states = automaton(patterns)
    certificate, _, prefactor = support(states, Q, 48)
    assert len(patterns) == 1207 and len(states) == 2899
    assert certificate == R and prefactor**2 < 2
    assert "all_odd_return_in 30 moves" in (ROOT/"output/research73/growth/i32.txt").read_text()
    audit = (OUT/"w32-audit.txt").read_text()
    root,maximum = map(int,re.search(r"root (\d+) largest (\d+)",audit).groups())
    assert X == Fraction(re.search(r"min (\d+/\d+)",audit)[1])
    factor = Fraction(root,maximum)/X**19
    assert factor >= Fraction(1,2**71)
    assert X>4 and X**18<2**71 and X<12
    assert (R*R/X)**16 < Q**9 and Q**3 > 4
    assert (R*R/OLD)**16 > Q**9
    data=dict(old_growth=str(OLD), new_growth=str(X), odd_return_sites=36,
              root=root, maximum_odd=maximum, uniform_power_of_two=-71,
              old_rational_exponent="4/7", new_rational_exponent="9/16",
              prescribed_error_exponent="3/2", Q=str(Q), r=str(R), examples=[])
    for n in (30,1000,20000):
        for eps in (Fraction(1,100),Fraction(1,10**6)):
            old=width(n,eps,OLD_FACTOR*OLD**n)
            new=width(n,eps,X**n/2**71)
            data["examples"].append(dict(order=n,epsilon=str(eps),old_sufficient_width=old,
                                         new_sufficient_width=new,combined_width=min(old,new)))
    return data


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--report",type=Path)
    parser.add_argument("--order",type=int)
    parser.add_argument("--epsilon",type=Fraction)
    parser.add_argument("--count",action="store_true",help="Also execute the exact confined transfer for the selected width")
    args=parser.parse_args()
    if (args.order is None) != (args.epsilon is None):
        parser.error("Supply both --order and --epsilon")
    if args.count and args.order is None:
        parser.error("--count requires --order and --epsilon")
    if args.order is not None:
        if args.order < 1 or not 0 < args.epsilon < 1 or args.report:
            parser.error("Use a positive order, accuracy in (0,1), and no --report")
        n, eps = args.order, args.epsilon
        denominator = max(OLD_FACTOR*OLD**n,X**n/2**71)
        chosen = width(n,eps,denominator)
        answer = dict(order=n, relative_error_at_most=str(eps), sufficient_width=chosen)
        if args.count:
            from strip_error import connected_strip
            answer["lower_count"] = connected_strip(chosen,(n,))[n]
        print(json.dumps(answer))
        sys.exit(0)
    data=check()
    result=json.dumps(data,indent=2)
    print(result)
    if args.report:
        args.report.write_text(result+"\n",encoding="utf-8")
    else:
        assert data == json.loads((OUT/"accuracy.json").read_bytes())
