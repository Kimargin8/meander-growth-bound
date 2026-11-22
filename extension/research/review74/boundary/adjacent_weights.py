"""Exact adjacent-gap contributions from outside walks of unrestricted length."""
from fractions import Fraction
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def weights(width, x):
    """x is the squared site parameter; gaps are numbered from 0."""
    if width < 2 or width % 2 or x <= 4:
        raise ValueError("Expected a positive even width and x greater than four")
    m = width + 3
    # det(sI-P_j) = s^(j mod 2) * reduced[j], with s^2 = x.
    reduced = [1, 1]
    for j in range(2, m+1):
        reduced.append((1 if j % 2 else x)*reduced[-1] - reduced[-2])
    return [Fraction(reduced[i]*reduced[m-i-2], reduced[m])
            for i in range(1, width+1)]


def check():
    # Solve for a source at each internal right endpoint independently.
    # At site parameter 4 the matrix is rational and strictly diagonally dominant.
    for w in range(2, 33, 2):
        m = w+3
        for i, expected in enumerate(weights(w,16), 1):
            diagonal = [Fraction(4)] * m
            rhs = [Fraction(j == i+1) for j in range(m)]
            for j in range(1,m):
                diagonal[j] -= 1/diagonal[j-1]
                rhs[j] += rhs[j-1]/diagonal[j-1]
            for j in reversed(range(m)):
                rhs[j] = (rhs[j] + (rhs[j+1] if j+1 < m else 0))/diagonal[j]
            assert rhs[i] == expected
    result = weights(32, 12)
    assert result == result[::-1] and min(result) > Fraction(1,9) > Fraction(7,64)
    return dict(width=32, squared_parameter=12,
                minimum=str(min(result)), maximum=str(max(result)),
                weights_by_internal_edge=list(map(str,result)),
                checked_edges=sum(range(2,33,2)))


if __name__ == "__main__":
    data = check()
    print(json.dumps(data,indent=2))
