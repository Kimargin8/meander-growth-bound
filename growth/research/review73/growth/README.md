# Root certificates for strip-confined meanders

The width is the total number of upper and lower road ends exposed by a vertical cut. Both banks and all river sites are labelled. An order-`n` meander has `2n` sites and one road component. Let `C_n(W)` count the meanders with width at most `W`.

The width-32 certificate places its growth constant between the exact fractions

\[
\frac{18325770610324}{1568251052756}
\ \le\ \mu_{32}\ \le\ 
\frac{15264493432827}{1305926315603}.
\]

In particular the unrestricted growth constant exceeds 11.685482. The upper endpoint is less than 11.688634 for this strip. An independently written generator reproduced all six exact endpoints from the three witness vectors.

## Why each verified fraction is a lower bound

After the first compulsory pair of openings, the frontier is a matching with two ends and the slit between them. Write `o` for this state. In the single-site transfer `T`, a state consists of a noncrossing matching of the exposed ends and the number of ends above the slit. Each permitted letter has weight one. The two shifts move the slit one position in opposite directions. A simultaneous opening inserts two paired ends at the slit. A simultaneous closure removes the two nearest ends and joins their former partners. It is permitted only when the two removed ends belonged to different paths. The width restriction applies to every single site. Returning to the root permits the final compulsory closure. Hence

\[
C_n(W)=(T^{2n-2})_{oo}.
\]

The bank-exchange involution carries a matching to its reverse reflection and sends a slit at `u` among `2k` ends to `2k-u`. Keeping representatives with `u≤k` gives an exact quotient for root counts. When two physical moves have the same quotient target their weights add. In particular the central slit still has two shift moves. Slit parity alternates at each step. Let `A=T²` on the odd slit class.

All physical states communicate within the width bound. Slit shifts can place the slit anywhere on the current matching. Every noncrossing matching has an adjacent paired ear. Removing such an ear in reverse time gives a predecessor under opening, and repeating this reaches the single pair. In the other direction, a matching with at least four ends has two cyclically consecutive, linearly adjacent ends which are not partners. A legal closure there reduces the width by two. The slit can be shifted to the root after the last reduction. The root has a one-slot return, so paths of one parity can be padded to a common length.

For a positive integer vector `v` the programs compare every quotient `(Av)_i/v_i` by integer cross-multiplication. If the minimum is `L` and the maximum is `U`, then `L v≤Av≤U v`. The positive-support test is essential for comparison with the root count. The folded programs propagate the root coordinate in Boolean arithmetic. For a common integer `r` they certify `(A^r e_o)_i≥1` on every odd state. With `K=max_i v_i` this gives

\[
(A^{m+r})_{oo}\ge (A^m v)_o/K \ge L^m v_o/K .
\]

The unrestricted growth constant is therefore at least `L`. The positive upper vector bounds the growth in this particular strip by `U`. The unfolded program checks the root's communicating class via two breadth-first searches. Its positive vector on both parities restricts to the odd class. Finite strong connectivity and the two root-returning shift pairs let one pad all returns from odd states to a common even length.

## Lower bounds at every finite order

For the width-32 witness after 300 rounds, every odd state can return to the root in 15 slots. The finite verifier computes the full vector `c=A^15 e_o` exactly. With `K=v_o min_i c_i/v_i`, the componentwise inequalities give `M_n≥K L^(n−16)` for `n≥16`. Recording the exact coordinates gives a stronger coefficient than `v_o/max v`.

To obtain a single bound of the form `M_n≥c L^n` for every `n≥1`, propagate up to 40 slots. When a layer grows large, divide its entries by a common power of two and round down. If the accumulated power is `E` and the integer layer is `x`, then `A^m e_o≥2^E x`. At 40 slots the tail coefficient is `2^E v_o(min_i x_i/v_i)/L^41`. For the remaining orders 1 through 40 the root entries give certified direct lower bounds. Taking the minimum of this finite list and the tail coefficient establishes the all-orders inequality. The receipt records exact fractions and short decimal lower bounds for both the full rational `L` and its two-decimal downward rounding.

## Witness generation and arithmetic

The initial odd witness is constant. Each iteration applies the two single-site transfers and scales the resulting integer vector. Rounding is used only to choose the *next* witness. Each printed pair of bounds is verified on the old witness using the unscaled integer image. The vector is reconstructed deterministically from the source and the stated number of iterations.

Before a transfer every witness entry is at most `2^52`. A single site has at most four moves, including multiplicity, so every measured numerator is at most `16·2^52`. Sums fit in 64 bits, while exact ratio comparisons and adaptive rescaling use 128-bit products. Slit positions on the width-32 folded graph fit in 32 bits. The largest matching rank stored in three bytes belongs to width 30 and is smaller than `C_15=9,694,845`. Input width and state-count bounds are checked during graph construction.

The two final candidates were fixed in advance as `2Av−23v` and `5Av−58v`. The original final lower ratio exceeds 11.6, so both candidates are positive. Their images follow from the same `Av` and `A²v`. The latter has entries at most `2^60`. Every signed subtraction is checked for positivity, and the largest cross-products fit within 128 bits. Both candidates are retained in the transcript. The first gives the tightest upper bound and the second the best lower bound.

`strip_bound.cpp` constructs the full labelled graph. `folded_bound.cpp` reconstructs the bank-exchange quotient from the established confined transfer and stores the outgoing edges. `implicit_bound.cpp` calculates those edges when needed and keeps both alternating slit parities in one value array. Its optional checkpoint stores only the odd entries. W14 checkpoints were restarted independently and reproduced the same final fractions.

## Reproduction

From the repository root, run the following commands in a 64-bit Linux environment with Python 3.10+, GCC with OpenMP, and GNU C++17 support for unsigned 128-bit integers.

```text
python3 research/review73/growth/verify.py
python3 research/review73/growth/verify.py --small-native --threads 2
python3 research/review73/growth/verify.py --replay w26 --replay f28 --threads 12
python3 research/review73/growth/verify.py --replay f30 --threads 12
python3 research/review73/growth/verify.py --replay i32 --threads 14
python3 research/review73/growth/verify.py --edges32 --threads 6
python3 research/review73/growth/verify.py --finite32 --threads 14
```

The default command checks the exact arithmetic in the saved transcripts, input hashes, and every literal meander of orders 1 through 8. The second command checks all three native methods against a literal pair-of-matchings enumeration for seven small widths. The replay commands rebuild the programs, regenerate the witnesses and repeat every componentwise comparison. The resulting rational pairs must match the stored transcripts exactly. The upper bounds in the transcripts concern individual strips. The unrestricted claim uses the lower bounds and the inclusion of each strip ensemble.

The all-edge audit decodes every resulting matching as an involution on labelled ends. Its two primary tests explicitly perform insertion of a pair and reconnection of the former partners after a legal closure. The finite check regenerates the last width-32 witness and its temporary 3,088,468,656-byte checkpoint, checks every odd coordinate in the root-dominated vector, and verifies all root layers up to slot 40.

The archive needs the three transfer C++ sources, the focused checks `check_strip_edges.cpp` and `finite_bound.cpp`, the independent `verify_vector.cpp` with its transcript `vector32-final.txt`, `verify.py`, `research/review57/strip/confined.cpp`, the project snapshot `references/oeis/2026-11-15/b005315.txt`, and the named transcript files with their JSON receipt under `output/research73/growth`. The build products and private regeneration checkpoints can be omitted.

The finite frontier transfer is an established enumeration method. Albert and Paterson proved an unrestricted lower bound of 11.380 in Theorem 1.1 of *Bounds for the growth rate of meander numbers*. The project compares the new integer certificates with that explicit published bound. Worldwide priority for the strongest explicit lower constant is a separate bibliographical question. [Jensen](https://doi.org/10.1088/0305-4470/33/34/301) · [Albert and Paterson](https://doi.org/10.1016/j.jcta.2005.02.006)
