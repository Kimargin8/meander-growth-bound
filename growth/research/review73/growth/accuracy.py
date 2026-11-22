"""Choose a physical strip attaining a prescribed relative error in ordinary M_n.

Run without arguments to verify the supporting data and sample width choices.
With --order N --epsilon P/Q, print the least sufficient even width for the
explicit bound below (full width is always permitted). Only standard Python.
"""
from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
STAR = ROOT / "research/review73/star"
OUT = ROOT / "output/research73/growth"
sys.path.insert(0, str(STAR))
from approximation import cycle_language
from height_envelope import automaton, support

Q = Fraction(159, 100)
RHO = Fraction(3898274927, 10**9)
LOWER = Fraction(583, 50)    # 11.66
PARITY_STEPS = 15
SCALE = 2**52


def finite_denominator():
    """Power-of-two prefactor valid at every order for base 11.68."""
    text = (OUT / "finite32.txt").read_text(encoding="ascii")
    base = Fraction(292, 25)
    candidates = {}
    levels = {int(m): (int(r) << int(e)) for m, e, r in
              re.findall(r"level (\d+) exponent (\d+) root (\d+) maximum", text)}
    assert set(levels) == set(range(41))
    exponent40 = int(re.search(r"level 40 exponent (\d+)", text)[1])
    for a, b, lower, root in re.findall(
            r"candidate (\d+) (\d+) lower (\d+/\d+) upper \d+/\d+ root (\d+)", text):
        assert Fraction(lower) > base
        factor = Fraction(re.search(rf"domination 40 {a} {b} (\d+/\d+)", text)[1])
        tail = (2**exponent40 * int(root)) * factor / base**41
        candidates[(a, b)] = min([tail]+[levels[m]/base**(m+1) for m in range(40)])
    assert len(candidates) == 3
    best = max(candidates, key=candidates.get)
    c = candidates[best]
    bits = max(0, c.denominator.bit_length()-c.numerator.bit_length())
    if (c.numerator << bits) < c.denominator:
        bits += 1
    assert Fraction(1, 2**bits) <= c and (bits == 0 or 2**(bits-1)*c < 1)
    return base, bits, best


def width(order, epsilon, uniform=None):
    """Least even W certified by the union bound, or full width 2n."""
    if order < 1 or not 0 < epsilon < 1:
        raise ValueError("require positive order and 0 < epsilon < 1")
    if order <= PARITY_STEPS:
        return 2*order
    # v(o)/max(v) >= 2^-52. C, the DFA prefactor, satisfies C^2 < 2.
    if uniform is None:
        denominator = LOWER**(order-1-PARITY_STEPS) / SCALE
    else:
        base, bits, _ = uniform
        denominator = base**order / 2**bits
    target = (2*order-1)*2 * RHO**(2*order) / (epsilon*denominator)
    left, right = -1, order
    while right-left > 1:
        middle = (left+right)//2
        if Q**(2*middle) >= target:
            right = middle
        else:
            left = middle
    if right < order:
        assert Q**(2*right) >= target
    if right > 0:
        assert Q**(2*(right-1)) < target
    return 2*right


def verify():
    text = (OUT / "i32.txt").read_text(encoding="ascii")
    assert "width 32 cells 810375650 root 1 odd 386058577 even 424317073" in text
    assert "all_odd_return_in 30 moves" in text
    tested = [(int(m), Fraction(v)) for m, v in
              re.findall(r"round (\d+) min (\d+/\d+) row", text)]
    # Use the first supporting row, so continuing the witness search does
    # not alter the result or the parameters of the accuracy theorem.
    iteration, witness = next((m, r) for m, r in tested if r > LOWER)
    words, counts = cycle_language(10)
    rows = automaton(words)
    rho, _, constant = support(rows, Q, 48)
    assert len(rows) == 2899 and sum(counts) == 1207 and rho == RHO
    assert constant**2 < 2
    assert (RHO**2/LOWER)**7 < Q**4
    assert Q**3 > 4  # log_Q(2) < 3/2
    uniform = finite_denominator()
    examples = [{"order": n, "epsilon": str(e), "elementary_width": width(n, e),
                 "certified_width": width(n, e, uniform)}
                for n in (15, 16, 30, 32, 39, 100, 1000, 10000)
                for e in (Fraction(1, 100), Fraction(1, 10**6))]
    return {"strip_width": 32, "comparison_lower": str(LOWER),
            "first_supporting_iteration": iteration, "integer_row_minimum": str(witness),
            "common_return_slots": PARITY_STEPS, "universal_integer_scale": SCALE,
            "uniform_finite_base": str(uniform[0]), "uniform_prefactor_power_of_two": -uniform[1],
            "uniform_witness": list(uniform[2]),
            "avoidance_states": len(rows), "forbidden_patterns": len(words),
            "Q": str(Q), "rho_upper": str(rho), "prefactor_squared_strict_upper": 2,
            "examples": examples}


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("--order", type=int)
    p.add_argument("--epsilon", type=Fraction)
    p.add_argument("--count", action="store_true", help="compute the constrained count as well as the sufficient width")
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    if (a.order is None) != (a.epsilon is None):
        p.error("--order and --epsilon go together")
    if a.count and a.order is None:
        p.error("--count requires --order and --epsilon")
    if a.order is not None:
        selected = width(a.order, a.epsilon, finite_denominator())
        if a.count:
            from strip_error import connected_strip
            answer = connected_strip(selected, (a.order,))[a.order]
            print(json.dumps({"order": a.order, "width": selected, "lower_count": answer,
                              "relative_error_at_most": str(a.epsilon)}))
        else:
            print(selected)
    else:
        result = verify()
        cert = OUT / "accuracy.json"
        if a.write:
            temp = cert.with_suffix(".pending")
            temp.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
            temp.replace(cert)
        else:
            assert result == json.loads(cert.read_bytes())
        print("Verified growth witness, avoidance bound, precision exponent and", len(result["examples"]), "finite width choices.")
