"""Width and bounded-strip approximation from the ten-letter cycle language.

All bank labels and the reading origin are fixed. A closed component on a
proper interval is incompatible with a connected full diagram. Both sides of
each cut therefore avoid every subword with a self-contained component. For the
full scan the same statement applies for n>5. The marking Q,1,1,1/Q gives
Q^((a+b)/2) on each side; reversing the suffix exchanges its two directions.

Let r be certified by Tv <= r v for a positive vector v. Dropping matching
requirements across the cut and using the same marking on both sides bounds
the nth-root limsup of the fraction with height at least c*n by
r^2/(mu*Q^c), where lim M_n^(1/n) >= mu. Union over all cuts and cyclic
origins costs only a polynomial. At all remaining cuts every side of length
ell >= n has d/ell > (1-c)/2 since 2d=ell-(a+b)+j and j >= 0.

At width W=2*floor(c*n/2) the exact confined transfer uses O*(2^(c*n))
bit time and space. Its relative error has the same exponential upper rate;
strip_error.py supplies a separate polynomial finite error bound at n>4
using factors of length at most eight. Neither calculation is an exact full
meander counter. The reported c is the best per-mille certificate on the
tested rational grid. No continuous global optimization is assumed.
"""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

from height_envelope import automaton, internal_cycle, rate_ceiling, support


SCENARIOS = (
    ("11.38", "Albert and Paterson, Theorem 1.1", "published lower bound"),
    ("11.44", "Width-26 lower certificate, growth/w26.txt, round 250", "project lower bound"),
    ("67010071789518/5852492708823", "Width-26 lower certificate, growth/w26.txt, round 250", "project lower bound"),
    ("11.50", "Width-28 lower certificate, growth/f28.txt", "project lower bound"),
    ("11.508", "Width-28 lower certificate, growth/f28.txt", "project lower bound"),
    *((s, "See the separate growth-rate certificate for applicability", "parametric") for s in
      ("11.53", "11.54", "11.57", "11.58", "11.60", "11.61", "11.65", "11.67", "11.70")),
)
SIMPLE_HEIGHTS = (Fraction(4,7), Fraction(58,100), Fraction(59,100), Fraction(3,5))


def cycle_language(limit):
    """A pair of dense byte arrays replaces the collection of forbidden words."""
    previous = bytearray(1)
    patterns, distribution = [], []
    for length in range(1, limit + 1):
        size = 4**length
        mask = size // 4 - 1
        current = bytearray(size)
        start = len(patterns)
        for value in range(size):
            if previous[value >> 2] or previous[value & mask]:
                current[value] = 1
            else:
                word = bytes(value >> shift & 3 for shift in range(2*length-2, -1, -2))
                if internal_cycle(word):
                    current[value] = 1
                    patterns.append(word)
        distribution.append(len(patterns) - start)
        previous = current
    known = set(patterns)
    assert all(bytes(3-x for x in reversed(w)) in known for w in patterns)
    assert all(bytes((x % 2)*2 + x//2 for x in w) in known for w in patterns)
    return patterns, distribution


def verify_transitions(rows, patterns, limit):
    """Check the DFA against a direct local graph test on every test word."""
    visited = 0
    stack = [(0, b"")]
    while stack:
        state, word = stack.pop()
        if word:
            assert (state < 0) == internal_cycle(word)
            visited += 1
        if len(word) < limit:
            for c in range(4):
                stack.append((rows[state][c] if state >= 0 else -1, word+bytes([c])))
    for word in patterns:
        assert internal_cycle(word)
        assert not internal_cycle(word[:-1]) and not internal_cycle(word[1:])
        state = 0
        for c in word[:-1]:
            state = rows[state][c]
            assert state >= 0
        assert rows[state][word[-1]] == -1
    return visited


def minimum_permille(rho, q, lower):
    # Binary search over integers; each comparison is between exact powers.
    target = (rho*rho/lower)**1000
    lo, hi = 0, 2000
    assert q**lo <= target < q**hi
    while hi - lo > 1:
        m = (lo + hi)//2
        if q**m > target:
            hi = m
        else:
            lo = m
    return Fraction(hi,1000)


def required_growth(grid, c):
    """The best upper estimate for rho^2 / Q^c on the tested grid."""
    def power(row):
        q, rho, _ = row
        return (rho*rho)**c.denominator/q**c.numerator
    choice = min(grid, key=power)
    target = power(choice)
    scale = 1000000
    lo, hi = 0, 20*scale
    while hi - lo > 1:
        mid = (hi + lo)//2
        if Fraction(mid,scale)**c.denominator >= target:
            hi = mid
        else:
            lo = mid
    upper = Fraction(hi,scale)
    assert Fraction(lo,scale)**c.denominator < target <= upper**c.denominator
    return choice, upper


def run():
    patterns, histogram = cycle_language(10)
    rows = automaton(patterns)
    assert histogram[:8] == [0, 1, 0, 2, 0, 14, 0, 116]
    word_checks = verify_transitions(rows, patterns, 10)
    rounds = 48
    grid = []
    for percent in range(150,181):
        q = Fraction(percent,100)
        rho, _, constant = support(rows, q, rounds)
        grid.append((q, rho, constant))
    basic = {c:required_growth(grid,c) for c in SIMPLE_HEIGHTS}
    scenarios = []
    for raw, source, status in SCENARIOS:
        mu = Fraction(raw)
        candidates = []
        for q, rho, constant in grid:
            c = minimum_permille(rho, q, mu)
            candidates.append((c,q,rho,constant))
        c = min(row[0] for row in candidates)
        # At the same c, maximize the strictly certified margin. All rates
        # have the same rational exponent and can be compared by integer powers.
        def objective(row):
            _, q, rho, _ = row
            return (rho*rho/mu)**c.denominator / q**c.numerator
        _, q, rho, constant = min((row for row in candidates if row[0] == c), key=objective)
        assert (rho*rho/mu)**c.denominator < q**c.numerator
        thresholds = []
        for target in sorted({c, *SIMPLE_HEIGHTS}):
            tq, trho, _ = basic[target][0] if target in basic else (q,rho,constant)
            value = rate_ceiling(trho,tq,target,mu)
            strict = (trho*trho/mu)**target.denominator < tq**target.numerator
            thresholds.append(dict(c=str(target), Q=str(tq), rho_upper=str(trho), margin_proves_decay=strict,
                                   rate_upper=str(value)))
        scenarios.append(dict(mu_lower=raw, status=status, source=source, Q=str(q),
                              rho_upper=str(rho), side_prefactor=str(constant),
                              c=str(c), excess_fraction=str((1-c)/2), thresholds=thresholds))
    script = Path(__file__).resolve()
    return dict(maximum_factor_length=10, factor_counts_by_length=histogram,
                words_examined=sum(4**s for s in range(1,11)), exhaustive_dfa_checks=word_checks,
                minimal_factors_checked=len(patterns), states=len(rows),
                alphabet=["UU","UD","DU","DD"], positive_vector_powers=rounds,
                factor_list_sha256=sha256(b"\n".join(patterns)).hexdigest(),
                grid=[dict(Q=str(q), rho_upper=str(rho)) for q,rho,_ in grid],
                sufficient_growth_thresholds=[dict(c=str(c), Q=str(choice[0]), rho_upper=str(choice[1]),
                                                    sufficient_if_mu_strictly_above=str(level))
                                              for c,(choice,level) in basic.items()],
                source_hashes={f.name:sha256(f.read_bytes()).hexdigest() for f in
                               (script, script.with_name("height_envelope.py"))},
                previous_length_8=dict(mu_lower="11.38", c="627/1000", rate_upper="999665/1000000"),
                scenarios=scenarios)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    d = run()
    if args.write:
        pending = args.write.with_suffix(".pending")
        pending.write_text(json.dumps(d, indent=2) + "\n",encoding="utf8")
        pending.replace(args.write)
    else:
        path = Path(__file__).resolve()
        cert = path.parent / "certificates/approximation.json"
        if not cert.is_file():
            cert = path.parents[3] / "output/research73/star/approximation.json"
        assert d == json.loads(cert.read_bytes())
    print(d["factor_counts_by_length"], d["states"])
    for row in d["scenarios"]:
        print(row["mu_lower"], row["Q"], row["rho_upper"], row["c"], row["excess_fraction"], row["thresholds"])
