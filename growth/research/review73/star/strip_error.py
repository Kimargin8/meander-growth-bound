"""Certified finite errors for the bounded-width enumeration of closed meanders.

Bank labels and the left end of the full Dyck word are fixed. At each cut the
width is the sum of the two open-arc counts. For n>4 a connected pair of Dyck
words of length 2n has no closed component in any interval of at most eight
vertices. If A(n,W) counts pairs satisfying this avoidance and the width bound,
then 0 <= M_n-C_n(W) <= A(n,2n)-A(n,W): every omitted meander is one of the
counted omitted words. The difference uses the same full-word automaton at each
width and enforces the empty final boundary on both banks.

For arbitrary n use W = 2 floor(627 n/2000). A suppressed width is at least
0.627n; the height-envelope theorem bounds the limsup nth root of the omitted
fraction by 0.999665. The exact restricted count is consequently a deterministic
relative approximation with exponentially small error as n grows. There are at
most sum_{2k<=W} (2k+1) Catalan(k) boundary states, counting the position of the
bank separator. Each has four transitions and O(n)-bit weights. This gives
O*(2^(0.627n)) bit time and space, with polynomial factors suppressed. The
asymptotic guarantee and the finite error bounds here are separate.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from pathlib import Path
import json

from height_envelope import automaton, dyck, forbidden_words


ORDERS = (30, 32, 39)
STRIPS = {
    20: ("output/research54/open/width20.txt",
         "8c6c0e23a4be67cc2165179755a78fb0f2a8ba0701ea7eb2ad9f7dbe827cf69f",
         (1236843595100226212566968598, 146904568491240703228063665980, 2806714409015470527492501800811643708)),
    24: ("output/research54/open/width24.txt",
         "df06cf44a1204c2ce0987781be9646dee47eb9429304e4cf15d0b1439bdd0271",
         (1276790853241587768819270968, 154240822354907504454705971832, 3234578560933647646617744371546454424)),
    28: ("output/research57/strip/pilot28.txt",
         "6dc186b68a5d0e6b21c5f17b5113fed7d317d91f343426b037de43150dc1d18d",
         (1279735608701986504826648142, 154932979355595311922036454298, 3312550539438300827339835757135053338)),
    30: ("output/research57/open/reverse30_51.txt",
         "caac91c660b00d388336a600a52fbb5e4fd58bb480e2741ae29520a107a021c9",
         (1279827515088389385323812816, 154961299926222703848278437100, 3318878242255407752020181251152722640)),
    32: ("output/research57/strip/reverse32_58.txt",
         "e131ba2d02a35693f127fa140df0bd2bd1356fbb4d087f74c1355b91e0d4ab82",
         (1279841059394216005714762516, 154966251760992175445253465146, 3320585060697892097353682548983016820)),
}
EXACT = {30: 1279842888607976889797535044, 32: 154967080111946539652848747626}
PREVIOUS = (3321085593816506271022366812311479295, 3321091581248739556828790840719943659)
NATIVE_WIDTH16 = (985420425595107241004955780, 109631330133350665277027672958, 1607855259195515751622595041906305752)
REFERENCES = {
    "output/research57/upper39/capacity10.json": "2a368d8e11a19dfa3f696de4a38f37343fb77fdf4f33149eb8f0fcf8f51b0cb7",
    "references/oeis/2026-11-15/b005315.txt": "f4a7b3ca6fc1ac097fa2f6aaef2d8988974984c97a52917a85735d15e048b3f9",
    "research/review57/strip/confined.cpp": "41a09776c24b00216632f0150fd5e4801b43b10328a4daee92abb87813fd5e3e",
}


def source_hash(file):
    return sha256(file.read_bytes()).hexdigest()


def audit_inputs(root):
    for width, (name, digest, values) in STRIPS.items():
        file = root / name
        assert source_hash(file) == digest
        lines = file.read_text().splitlines()
        assert lines[0].startswith(f"# width {width} ")
        counts = {int(row[1]): int(row[2]) for line in lines if (row := line.split()) and row[0] == "C"}
        assert tuple(counts[n] for n in ORDERS) == values
    for name, digest in REFERENCES.items():
        assert source_hash(root / name) == digest
    saved = json.loads((root / "output/research57/upper39/capacity10.json").read_bytes())
    assert (saved["lower"], saved["upper"]) == PREVIOUS
    table = (root / "references/oeis/2026-11-15/b005315.txt").read_text().splitlines()
    terms = {int(row[0]): int(row[1]) for line in table if (row := line.split()) and row[0].isdigit()}
    assert all(terms[n] == value for n, value in EXACT.items())


def avoided_counts(rows, orders, width):
    """Count words exactly. All totals use Python integers."""
    target = 2 * max(orders)
    moves = [tuple((1 if not c & 2 else -1, 1 if not c & 1 else -1, j)
                   for c, j in enumerate(row) if j >= 0) for row in rows]
    layer = {(0, 0, 0): 1}
    totals = {}
    for step in range(1, target + 1):
        remaining = target - step
        following = {}
        for (a, b, state), count in layer.items():
            for da, db, dest in moves[state]:
                h, k = a + da, b + db
                if min(h, k) < 0 or max(h, k) > remaining or h + k > width:
                    continue
                key = h, k, dest
                following[key] = following.get(key, 0) + count
        layer = following
        if step % 2 == 0 and step // 2 in orders:
            totals[step // 2] = sum(v for (a, b, _), v in layer.items() if a == b == 0)
    return totals


def connected_strip(width, orders):
    """Exact physical transfer for small strip widths, including the width-16 input."""
    @cache
    def extensions(pairs, a):
        # Top stubs run outside-in, bottom stubs inside-out. The slit is at a.
        result = []
        s = len(pairs)
        if s < width:
            shifted = tuple(v + 2 if v >= a else v for v in pairs)
            result.append((shifted[:a] + (a + 1, a) + shifted[a:], a + 1))
        if a:
            result.append((pairs, a - 1))
        if a < s:
            result.append((pairs, a + 1))
        if 0 < a < s and pairs[a-1] != a:
            # DD reconnects the other ends of the arcs meeting the slit.
            i, j = pairs[a-1], pairs[a]
            joined = list(pairs)
            joined[i], joined[j] = j, i
            raw = joined[:a-1] + joined[a+1:]
            result.append((tuple(v-2 if v > a else v for v in raw), a-1))
        return tuple(result)

    layer = {((), 0): 1}
    counts = {}
    for step in range(1, 2 * max(orders) + 1):
        if step % 2 == 0 and step // 2 in orders:
            # Only a final DD from this two-stub state closes the sole component.
            counts[step // 2] = layer.get(((1, 0), 1), 0)
        following = {}
        for (pairs, a), count in layer.items():
            for state in extensions(pairs, a):
                following[state] = following.get(state, 0) + count
        layer = following
    return counts


def literal_checks(rows, patterns):
    """Enumerate actual closed diagrams and words by their maximum width."""
    bad = set(map(tuple, patterns))
    sizes = sorted(set(map(len, bad)))
    checked = []
    for n in (5, 6):
        def arcs(w):
            partners = [0] * len(w)
            stack = []
            for i, bit in enumerate(w):
                if bit == 0:
                    stack.append(i)
                else:
                    j = stack.pop()
                    partners[i], partners[j] = j, i
            return partners
        words = [(w, arcs(w)) for w in dyck(n)]
        connected, allowed = Counter(), Counter()
        for u, p in words:
            for v, q in words:
                word = tuple(2 * x + y for x, y in zip(u, v))
                h = maximum = 0
                for letter in word:
                    h += 2 * ((letter == 0) - (letter == 3))
                    maximum = max(maximum, h)
                admissible = not any(word[i:i+k] in bad for k in sizes for i in range(len(word)-k+1))
                if admissible:
                    allowed[maximum] += 1
                i, visited = 0, 0
                while True:
                    i = q[p[i]]
                    visited += 2
                    if i == 0:
                        break
                if visited == 2*n:
                    assert admissible
                    connected[maximum] += 1
        assert sum(connected.values()) == {5: 262, 6: 1828}[n]
        for width in range(0, 2*n+1, 2):
            expect_a = sum(v for h, v in allowed.items() if h <= width)
            expect_c = sum(v for h, v in connected.items() if h <= width)
            assert avoided_counts(rows, (n,), width)[n] == expect_a
            assert connected_strip(width, (n,))[n] == expect_c
        checked.append(dict(n=n, dyck_pairs=len(words)**2, connected_by_max_width=sorted(connected.items()),
                            allowed_by_max_width=sorted(allowed.items())))
    return checked


def run():
    patterns, histogram = forbidden_words(8)
    rows = automaton(patterns)
    small = literal_checks(rows, patterns)
    full = avoided_counts(rows, ORDERS, 2*max(ORDERS))
    confined = {16: connected_strip(16, ORDERS)}
    assert tuple(confined[16][n] for n in ORDERS) == NATIVE_WIDTH16
    confined.update({w: dict(zip(ORDERS, values)) for w, (_, _, values) in STRIPS.items()})
    intervals = []
    best = dict.fromkeys(ORDERS)
    last_confined = dict.fromkeys(ORDERS, 0)
    last_allowed = dict.fromkeys(ORDERS, 0)
    for width, counts in sorted(confined.items()):
        partial = avoided_counts(rows, ORDERS, width)
        for n in ORDERS:
            lo = counts[n]
            tail = full[n] - partial[n]
            upper = lo + tail
            assert last_confined[n] <= lo <= partial[n] <= full[n] and last_allowed[n] <= partial[n]
            assert lo - last_confined[n] <= partial[n] - last_allowed[n]
            last_confined[n], last_allowed[n] = lo, partial[n]
            row = dict(n=n, width=width, confined_meanders=lo, confined_avoiding_pairs=partial[n],
                       all_avoiding_pairs=full[n], tail_upper=tail, interval=[lo, upper],
                       relative_error_upper=str(Fraction(tail, upper)))
            if n in EXACT:
                assert lo <= EXACT[n] <= upper
                row["true_tail"] = EXACT[n] - lo
            else:
                assert lo <= PREVIOUS[1] and upper >= PREVIOUS[0]
            intervals.append(row)
            if best[n] is None or upper < best[n]["interval"][1]:
                best[n] = row
    lo39 = max(PREVIOUS[0], max(row["interval"][0] for row in intervals if row["n"] == 39))
    upper39 = best[39]["interval"][1]
    assert lo39 <= min(upper39, PREVIOUS[1])
    this = Path(__file__).resolve()
    return dict(convention="Fixed bank labels and scan origin; maximum over all cuts of a+b; ordinary closed connected meanders.",
                avoidance=dict(maximum_factor_length=8, minimal_factor_counts=histogram, states=len(rows)),
                exact_reference=EXACT, previous_order39_interval=PREVIOUS,
                source_hashes={p.name:source_hash(p) for p in (this, this.with_name("height_envelope.py"))},
                prior_full_file_hashes={**{name:h for name,h,_ in STRIPS.values()}, **REFERENCES},
                width16_independent_check=dict(source="research/review57/strip/confined.cpp",
                                               arguments=[16,39],orders=ORDERS,counts=NATIVE_WIDTH16),
                asymptotic=dict(width="2 floor(627 n/2000)", relative_error_nth_root_limsup_upper="999665/1000000",
                                bit_time_and_space="O*(2^(0.627 n))", source="height_envelope.py"),
                small_checks=small, intervals=intervals,
                order39_comparison=dict(best_width=best[39]["width"], best_upper=upper39,
                                        improves_previous_upper=upper39 < PREVIOUS[1],
                                        interval_intersection=[lo39, min(upper39, PREVIOUS[1])]))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path)
    parser.add_argument("--audit-root", type=Path, help="Check the original strip and reference tables by full-file hash and selected terms.")
    opts = parser.parse_args()
    if opts.audit_root:
        audit_inputs(opts.audit_root)
    data = run()
    if opts.write:
        pending = opts.write.with_suffix(".pending")
        pending.write_text(json.dumps(data, indent=2) + "\n", encoding="utf8")
        pending.replace(opts.write)
    else:
        script = Path(__file__).resolve()
        file = script.parent / "certificates/strip-error.json"
        if not file.is_file():
            file = script.parents[3] / "output/research73/star/strip-error.json"
        assert json.loads(file.read_bytes()) == json.loads(json.dumps(data))
    print(json.dumps(data["order39_comparison"]))
