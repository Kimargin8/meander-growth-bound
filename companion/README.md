# Shortest and near-shortest halves of closed plane meanders

Iaroslav Argunov. The [paper](shortest_meander_halves.pdf) gives exact multiplicities for shortest cut types, a counting algorithm at each fixed excess, and height bounds for a uniform closed meander. [LaTeX source](shortest_meander_halves.tex).

The top-level [verification program](../verify.py) runs these checks in a fresh Python process for each file. A check of one statement can also be run separately from this directory.

| Statement or computation | Program |
|---|---|
| Multiplicities and generating functions by cut type | [minimum_caps.py](minimum_caps.py) |
| Realization of a prescribed type at fixed excess | [check_excess_caps.py](check_excess_caps.py) and [excess_caps.py](excess_caps.py) |
| Sublinear excess | [check_excess_growth.py](check_excess_growth.py) |
| Probability at the central cut | [check_central_excess.py](supplement/research/review73/star/check_central_excess.py) |
| Simultaneous bounds for all cuts | [height_envelope.py](height_envelope.py) |
| Error estimates at orders 30, 32 and 39 | [cut_excess_bound.py](cut_excess_bound.py) and [strip_error.py](strip_error.py) |
| Height bound using cycle factors through length ten | [approximation.py](approximation.py) |
| Gluing at all small cuts and the obstruction for strictly shorter caps | [all_boundaries.py](all_boundaries.py) and [strict_order_kernel.py](supplement/research/review73/wholecut/strict_order_kernel.py) |

Ordinary Python 3.10 or later suffices. The programs use exact integers or rational comparisons. The recorded growth-rate input used in the height argument is available in the [main paper's certificates](../growth/output/research73/growth/growth-bounds.json).

The finite error check also reads the saved transfer counts for widths 20, 24, 28, 30 and 32 and compares their byte digests with its reference list. From this directory run `python3 -B strip_error.py --audit-root ../growth`. The [integer transfer](../growth/research/review57/strip/confined.cpp) accepts a width and a maximum order, so each confined count in the table can be regenerated.
