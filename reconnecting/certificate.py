"""Lightweight independent checks for the reconnecting-returns supplement."""
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import subprocess
import sys
from kernel import tables, literal
from upper_weights import weights
from words import computed, direct

HERE = Path(__file__).resolve().parent
DATA = HERE / "certificates"


def fixed_walk(width, start, target):
    counts = [0] * (width + 1)
    counts[start] = 1
    value = Q(0)
    for length in range(7):
        if length % 2 == 0:
            value += Q(counts[target], 12 ** (length // 2))
        counts = [a + b for a, b in zip([0] + counts[:-1], counts[1:] + [0])]
    return value


def main():
    if not __debug__:
        raise RuntimeError("Run in ordinary Python with assertions enabled")
    x = Q(1461549008357498, 124174541777295)
    u = Q(3431031596329, 292144566372)
    assert Q("11.7701179922") < x < 12
    assert Q("11.73") < u < Q("11.7442937206")
    assert x - u > Q("0.02582427")
    # Bound the infinite-line unit recursion by rearrangement and squaring.
    t = Q(7, 6)
    assert Q(117, 10) <= (1+2*t*(t-1))**2 * (Q(117, 10)-4*t*t)
    lo, alt = tables(32)
    assert lo == alt == {1:[582,593]+[595]*28+[593,582],
                          2:[63]+[64]*29+[63], 3:[83]+[84]*28+[83]}
    finite = literal(32, 6, 83)
    for d, limit in [(1,47),(2,28),(3,43)]:
        for k, row in zip(lo[d], finite[d]):
            assert sum(Q(v, 12**(j//2)) for j,v in enumerate(row[:limit+1])) >= Q(k,4096)
    for gap in range(31):
        fixed = fixed_walk(34,gap+1,gap+1)*fixed_walk(36,gap+2,gap+4)*fixed_walk(34,gap+3,gap+3)/12**2
        assert fixed >= Q(5,4096)
    assert weights("continuants") == weights("walks") == [611,624]+[626]*28+[624,611]
    for rec in json.loads((DATA/"small-words.json").read_bytes()):
        counts = computed(rec["width"], rec["depth"], 7)
        assert counts == rec["meanders"] == [direct(rec["width"],rec["depth"],n) for n in range(1,8)]
    for typ, value in [("mixed",x),("unit",u)]:
        receipt = json.loads((DATA/(typ+".json")).read_bytes())
        reference = (DATA/(typ+"-audit.txt")).read_text()
        repeat = (DATA/(typ+"-rescan.txt")).read_text()
        for log in (reference, repeat):
            bound = re.search(r"min (\d+/\d+) max (\d+/\d+)",log)
            assert (Q(bound[1]),Q(bound[2]))==(Q(receipt["minimum"]),Q(receipt["maximum"]))
            assert "odd 386058577 even 424317073" in log
        assert Q(receipt["certified"])==value
    for q,r,power,num in [(Q(159,100),Q(19392772299,5000000000),17,9),
                           (Q(36,25),Q(18853838429,5000000000),27,14)]:
        assert (r*r/Q("11.7701"))**power < q**num
    assert Q(159,100)**3 > 4 and Q(36,25)**2 > 2
    print("Verified the coefficients, all-depth upper envelope, 63 small-word counts, two exact extrema, and approximation exponents.")


if __name__ == "__main__":
    main()
