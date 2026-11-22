"""Independently test the projected-prefix lemma and its rational prefactor."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path


def move(pair, gap, symbol):
    if symbol == "L":
        return (pair, gap-1) if gap else None
    if symbol == "R":
        return (pair, gap+1) if gap < len(pair) else None
    if symbol == "O":
        relabel = tuple(x+2 if x >= gap else x for x in pair)
        return relabel[:gap]+(gap+1,gap)+relabel[gap:], gap+1
    if not 0 < gap < len(pair) or pair[gap-1] == gap:
        return None
    left, right = pair[gap-1], pair[gap]
    joined = list(pair)
    joined[left], joined[right] = right, left
    return tuple(x-2 if x > gap else x for x in joined[:gap-1]+joined[gap+1:]), gap-1


def annotations(width, through):
    # Each entry is a labelled expanded walk from the root; multiple private
    # annotations of a single word remain separate entries in this list.
    layer = [("", (1,0), 1, None)]
    checked = ambiguous = odd_outside = 0
    for steps in range(through+1):
        multiplicities = Counter(word for word, _, _, _ in layer)
        for word, pairing, gap, target in layer:
            k = multiplicities[word]
            assert k <= (1 if target is None else 2)
            if target is not None and steps % 2 == 0:
                assert gap % 2 and len(pairing) == width+2
                if move(pairing,gap,"C") is None:
                    neighbour = gap+1 if gap < len(pairing)-1 else gap-1
                    assert move(pairing,neighbour,"C") is not None
                odd_outside += 1
        checked += len(multiplicities)
        ambiguous += sum(count == 2 for count in multiplicities.values())
        if steps == through:
            break
        following = []
        for word, pair, gap, target in layer:
            symbols = "LRCO" if target is None else ("LRC" if gap==target else "LR")
            for letter in symbols:
                state = move(pair,gap,letter)
                if state is None:
                    continue
                if target is None and letter == "O" and len(pair) == width:
                    for d in (-1,1):
                        if 0<=gap+d<=width:
                            following.append((word+letter,*state,gap+1+d))
                else:
                    following.append((word+letter,*state,None if letter=="C" else target))
        layer = following
    assert ambiguous and odd_outside
    return dict(width=width, steps=through, distinct_prefixes_checked=checked,
                ambiguous_prefixes=ambiguous, annotated_odd_outside_ends=odd_outside)


def constants(x, root, maximum):
    if not (root>0 and maximum>=root and x>4):
        raise ValueError("Positive coordinates and a squared site parameter greater than 4 required")
    assert x >= 8 and x*(x-4) >= 64
    factor = Fraction(root,maximum)/x**17
    p = 0
    while (2**p)*factor < 1 or 2**p < x**16:
        p += 1
    q, r = Fraction(159,100), Fraction(3898274927,10**9)
    crosses = (r*r/x)**25 < q**14
    assert q**3 > 4
    return dict(squared_parameter=str(x), root=root, maximum_odd=maximum,
                all_order_power_of_two=-p, projected_completion=32,
                q=str(q), r=str(r), exponent_14_over_25_supported=crosses)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--x",type=Fraction,required=True)
    parser.add_argument("--root",type=int,required=True)
    parser.add_argument("--maximum",type=int,required=True)
    parser.add_argument("--report",type=Path)
    args=parser.parse_args()
    data=constants(args.x,args.root,args.maximum)
    data["enumeration"]=[annotations(w,8) for w in (2,4,6)]
    result=json.dumps(data,indent=2)
    print(result)
    if args.report:
        args.report.write_text(result+"\n",encoding="utf-8")
