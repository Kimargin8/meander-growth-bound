"""Check exact elimination of a full outer strip layer at a rational parameter."""
from argparse import ArgumentParser
from hashlib import sha256
from math import comb
from pathlib import Path
from time import perf_counter
import json
import platform

ROOT = Path(__file__).resolve().parents[3]


def words(length, position=0, height=0, word=0):
    if position == length:
        yield word
        return
    if height < length - position:
        yield from words(length, position+1, height+1, word | (1 << position))
    if height:
        yield from words(length, position+1, height-1, word)


def partners(word, length):
    mate = [-1] * length
    stack = []
    for i in range(length):
        if word >> i & 1:
            stack.append(i)
        else:
            j = stack.pop()
            mate[i], mate[j] = j, i
    return mate


def remove(word, left):
    return (word & ((1 << left)-1)) | ((word >> (left+2)) << left)


def after_close(word, gap, mate):
    a, b = mate[gap-1:gap+1]
    assert a != gap
    joined = (word | (1 << min(a, b))) & ~(1 << max(a, b))
    return remove(joined, gap-1)


def continuants(size, p, q):
    result = [1, p]
    for _ in range(2, size+1):
        result.append(p*result[-1] - q*q*result[-2])
    return result


def solve_path(b, p, q, det):
    """Return det[m] * (p I - q P)^(-1) q b, using integer arithmetic."""
    forward = []
    previous = 0
    for i, item in enumerate(b):
        previous = q * (item*det[i] + previous)
        forward.append(previous)
    answer = [0] * len(b)
    right = 0
    for i in reversed(range(len(b))):
        value = det[-1]*forward[i] + q*det[i]*right
        right, remainder = divmod(value, det[i+1])
        assert remainder == 0
        answer[i] = right
    return answer


def trial_value(word, split, width):
    # Positive, non-symmetric test vector, determined without an exterior array.
    return 3 + ((word*227 + split*127 + width*773) % 997)


def scan(width, p, q):
    start = perf_counter()
    length = width + 2
    size = length + 1
    det = continuants(size, p, q)
    outputs = {}
    matchings = 0
    for word in words(length):
        mate = partners(word, length)
        b = [0] * size
        for i in range(1, length):
            if mate[i-1] != i:
                b[i] = trial_value(after_close(word, i, mate), i-1, width)
        y = solve_path(b, p, q, det)
        for i, yi in enumerate(y):
            neighbors = (y[i-1] if i else 0) + (y[i+1] if i < length else 0)
            assert p*yi-q*neighbors == q*det[-1]*b[i]
        for i in range(1, length):
            if mate[i-1] == i:
                key = (remove(word, i-1), i-1)
                assert key not in outputs
                outputs[key] = y[i]
        matchings += 1
    elapsed = perf_counter()-start
    assert matchings == comb(length, length//2)//(length//2+1)
    assert len(outputs) == (width+1)*comb(width,width//2)//(width//2+1)
    # Independent row-oriented check: physically insert the pair and address the
    # path Green function by its cofactor formula; do not use the forward solve.
    for word in words(width):
        mate_small = partners(word, width)
        for gap in range(width+1):
            move = lambda k: k if k<gap else k+2
            original = [move(k) for k in mate_small]
            large = original[:gap] + [gap+1,gap] + original[gap:]
            i = gap+1
            result = 0
            for j in range(1, length):
                a, b = large[j-1:j+1]
                if a == j:
                    continue
                remaining = [k for k in range(length) if k not in (j-1,j)]
                flags = 0
                for r, pos in enumerate(remaining):
                    partner = b if pos==a else a if pos==b else large[pos]
                    if partner>pos:
                        flags |= 1 << r
                cofactor = q**(1+abs(i-j))*det[min(i,j)]*det[size-1-max(i,j)]
                result += cofactor*trial_value(flags, j-1, width)
            assert outputs[word,gap] == result
    digest = sha256()
    for (word, gap), val in sorted(outputs.items()):
        digest.update(f"{word} {gap} {val}\n".encode("ascii"))
    return dict(width=width, outer_matchings=matchings, exterior_cells=matchings*size,
                returned_inner_rows=len(outputs), represented_durations="all",
                denominator=det[-1], numerators_sha256=digest.hexdigest(),
                streaming_apply_and_residual_check_seconds=round(elapsed,6))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    data = dict(parameter="7/2 per site", python=platform.python_version(),
                rows=[scan(w,7,2) for w in (2,4,6,8,10,12,14)])
    if args.report:
        args.report.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(data))


if __name__ == "__main__":
    main()
