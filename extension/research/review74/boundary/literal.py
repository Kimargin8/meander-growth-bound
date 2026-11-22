"""Compare excursion-polynomial counts with complete pairs of noncrossing matchings."""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "research/review74/excursions"))
from scan import matrices


@lru_cache(None)
def arches(n):
    if not n:
        return [()]
    result = []
    for k in range(n):
        right = 2*k+1
        for a in arches(k):
            for b in arches(n-k-1):
                result.append((right,)+tuple(i+1 for i in a)+(0,)+tuple(i+right+1 for i in b))
    return result


def interior_widths(pair):
    h = 0
    result = [h]
    for i, j in enumerate(pair):
        h += 1 if i < j else -1
        result.append(h)
    return result


def connected(a, b):
    pos, pairs = b[a[0]], 1
    while pos:
        pos = b[a[pos]]
        pairs += 1
    return 2*pairs == len(a)


def census(n, widths):
    pairs = arches(n)
    heights = [interior_widths(a) for a in pairs]
    counts = Counter()
    for i, a in enumerate(pairs):
        for j in range(i+1):
            if not connected(a, pairs[j]):
                continue
            levels = [x+y for x,y in zip(heights[i],heights[j])]
            peak = max(levels)
            for w in widths:
                if peak > w+2:
                    continue
                longest = 0
                run = 0
                for h in levels:
                    if h>w:
                        run += 1
                    elif run:
                        longest = max(longest, run+1)
                        run = 0
                for cutoff in (1,3,4,5,6):
                    if longest <= cutoff:
                        counts[w,cutoff] += 1 + int(i!=j)
    return counts


def verify(through=7):
    import numpy as np
    widths = (2,4,6,8,10)
    direct = {n: census(n,widths) for n in range(1,through+1)}
    rows = []
    for w in widths:
        edges = matrices(w,6)
        # scan.matrices lists (the sole size-two pairing, gap 0/1/2) first.
        root = 1
        for cutoff in (1,3,4,5,6):
            v = np.zeros(edges[1].shape[0],dtype=np.int64)
            v[root] = 1
            history = [v]
            for t in range(1,2*through-1):
                nxt = np.zeros_like(v)
                for d,m in edges.items():
                    if d<=min(cutoff,t):
                        nxt += m @ history[t-d]
                history.append(nxt)
            checks = [int(history[2*n-2][root]) for n in range(1,through+1)]
            expected = [direct[n][w,cutoff] for n in range(1,through+1)]
            assert checks == expected, (w,cutoff,checks,expected)
            rows.append(dict(width=w,cutoff=cutoff,counts=checks))
    return dict(orders=[1,through], comparisons=len(rows)*through, rows=rows)


if __name__ == "__main__":
    data=verify()
    (ROOT/"output/research74/boundary/literal.json").write_text(json.dumps(data,indent=2)+"\n")
    print("Agreement on",data["comparisons"],"complete enumeration counts.")
