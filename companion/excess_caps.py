"""Count loop-free physical halves of a prescribed boundary and length L+2*d.

The cut labels run from the outermost upper stub to the outermost lower stub.
No allocations depend on the number of shuffles of shortest components.
"""
from collections import defaultdict
from functools import cache
from itertools import combinations, permutations
from math import comb, prod
from operator import index


def positive_parts(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in positive_parts(total - first, parts - 1):
            yield (first, *rest)


def ordinary_graph(a, pairing, special):
    """The shortest ordinary half is a shuffle of whole Dyck excursions between through pairs."""
    size = len(pairing)
    top = [x for x in reversed(range(a)) if x not in special]
    bottom = [x for x in range(a, size) if x not in special]
    banks = (top, bottom)
    positions = [{x: i for i, x in enumerate(bank)} for bank in banks]
    through = [set(x for x in bank if (pairing[x] < a) != (x < a)) for bank in banks]

    def inside(bank, side):
        height = 0
        flags = [False]
        for i, x in enumerate(bank):
            if x not in through[side]:
                height += 1 if positions[side][pairing[x]] > i else -1
            flags.append(bool(height))
        return flags

    open_top, open_bottom = (inside(bank, side) for side, bank in enumerate(banks))
    # Each ordinary point is named by its stub, choosing the upper stub on a through path.
    def vertex(x):
        return pairing[x] if x >= a and pairing[x] < a else x

    def moves(state):
        i, h = state
        if i < len(top) and top[i] not in through[0] and not open_bottom[h]:
            yield (i + 1, h), vertex(top[i])
        if h < len(bottom) and bottom[h] not in through[1] and not open_top[i]:
            yield (i, h + 1), vertex(bottom[h])
        if i < len(top) and h < len(bottom) and pairing[top[i]] == bottom[h]:
            yield (i + 1, h + 1), top[i]

    # Use disjoint names for exposed stubs and for the new points on the line.
    # In each bank, two marked endpoints can be joined precisely when their masks
    # of open ordinary arcs agree. A set bit identifies an individual arc.
    masks_at_old = {}
    toggles = defaultdict(lambda: [0, 0])
    initial = []
    for side, bank in enumerate(banks):
        incident = {}
        for bit, stub in enumerate(bank):
            incident[stub] = 1 << bit
            toggles[vertex(stub)][side] ^= 1 << bit
        next_bit = len(bank)
        other_bank = banks[1-side]
        for stub in other_bank:
            mate = pairing[stub]
            if mate > stub and mate in positions[1-side]:
                bit = 1 << next_bit
                next_bit += 1
                toggles[vertex(stub)][side] ^= bit
                toggles[vertex(mate)][side] ^= bit
        mask = 0
        opening_order = range(a) if side == 0 else reversed(range(a, size))
        for stub in opening_order:
            if stub in special:
                masks_at_old[stub] = mask
            else:
                mask ^= incident[stub]
        initial.append(mask)
    start = (0, 0)
    end = (len(top), len(bottom))
    masks = {start: tuple(initial)}
    stack = [start]
    successors = {}
    while stack:
        state = stack.pop()
        successors[state] = []
        for nxt, point in moves(state):
            new = tuple(m ^ t for m, t in zip(masks[state], toggles[point]))
            successors[state].append(nxt)
            if nxt in masks:
                assert masks[nxt] == new
            else:
                masks[nxt] = new
                stack.append(nxt)
    assert masks[end] == (0, 0)
    # Use a topological traversal, including the empty interval from a state to itself.
    states = tuple(sorted(masks, key=lambda s: (sum(s), s)))

    @cache
    def path_counts(origin):
        counts = {origin: 1}
        for s in states:
            if s in counts:
                for nxt in successors[s]:
                    counts[nxt] = counts.get(nxt, 0) + counts[s]
        return counts
    return start, end, masks, masks_at_old, path_counts


def shape_orders(a, size, paths):
    """Choose the line order of the vertices of each oriented exceptional path."""
    # Old endpoint name -1-x denotes stub x. All new names are nonnegative.
    edges = [[], []]
    total = 0
    for x, y, length in paths:
        vertices = [-1-x, *range(total, total+length), -1-y]
        bank = int(x >= a)
        for u, v in zip(vertices, vertices[1:]):
            edges[bank].append((u, v))
            bank = 1-bank
        assert 1-bank == int(y >= a)
        total += length
    old_key = {-1-x: x-size for x in range(a)}
    old_key.update({-1-x: -1-x for x in range(a, size)})
    for order in permutations(range(total)):
        keys = {**old_key, **{x:i for i, x in enumerate(order)}}
        requirements = [[] for x in order]
        for bank, bank_edges in enumerate(edges):
            prior = []
            for u, v in bank_edges:
                if keys[u] > keys[v]:
                    u, v = v, u
                left, right = keys[u], keys[v]
                if any(i < left < j < right or left < i < right < j for i, j in prior):
                    break
                prior.append((left, right))
                # Both endpoints of an internal edge are markers. Relations to
                # future markers are tested when the later endpoint is visited.
                requirements[right].append((bank, u if u < 0 else keys[u]))
            else:
                continue
            break
        else:
            yield requirements


def count_near_minimum(a, pairing, d):
    """Count exactly one pairing with d pairs of vertices beyond its shortest length."""
    a, d = index(a), index(d)
    p = tuple(map(index, pairing))
    size = len(p)
    if not 0 <= a <= size or d < 0 or size % 2:
        raise ValueError("Expected a nonnegative excess and an even number of exposed ends")
    if any(y < 0 or y >= size or x == y or p[y] != x for x, y in enumerate(p)):
        raise ValueError("The boundary must be a perfect matching")
    stack = []
    for x, y in enumerate(p):
        if x < y:
            stack.append(y)
        elif stack.pop() != x:
            raise ValueError("The boundary must be noncrossing")
    if not size:
        return int(d == 0)  # formal empty half; no positive-length loop-free half
    pairs = [(x, y, 1 if x < a <= y else 2) for x, y in enumerate(p) if x < y]
    if d == 0:
        start, end, _, _, count = ordinary_graph(a, p, set())
        return count(start)[end]
    answer = 0
    for s in range(1, min(len(pairs), d)+1):
        for chosen in combinations(pairs, s):
            special = {v for x, y, _ in chosen for v in (x, y)}
            start, end, masks, old_masks, count = ordinary_graph(a, p, special)
            by_mask = [defaultdict(set), defaultdict(set)]
            for state, pair in masks.items():
                for side in (0, 1):
                    by_mask[side][pair[side]].add(state)
            for excesses in positive_parts(d, s):
                paths = [(x, y, L+2*k) for (x, y, L), k in zip(chosen, excesses)]
                for requirements in shape_orders(a, size, paths):
                    history = []
                    def interleave(state, weight):
                        k = len(history)
                        if k == len(requirements):
                            return weight * count(state).get(end, 0)
                        reachable = count(state)
                        candidates = reachable.keys()
                        for bank, endpoint in requirements[k]:
                            mask = old_masks[-1-endpoint] if endpoint < 0 else masks[history[endpoint]][bank]
                            candidates = candidates & by_mask[bank][mask]
                        subtotal = 0
                        for nxt in candidates:
                            history.append(nxt)
                            subtotal += interleave(nxt, weight * reachable[nxt])
                            history.pop()
                        return subtotal
                    answer += interleave(start, 1)
    return answer

