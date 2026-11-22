# A certified return of one fixed state

At width 26 a single return calculation proves **\(\mu>11.421465\)**. The transfer contains 14,061,140 even-slit and 13,027,729 odd-slit physical states. Independent enumeration of full diagrams reproduces all 42 width-restricted counts at widths 2 through 12 and orders 1 through 7.

| Two-site slots \(m\) | Accumulated binary exponent \(E_m\) | Root mantissa \(x_m(o)\) | Strict lower bound on \(\mu\) |
|---:|---:|---:|---:|
| 1,400 | 4,914 | 10,407,548,347,941 | 11.382493 |
| 3,100 | 10,894 | 11,568,213,483,477 | 11.421465 |

For the last row, the integer inequality is \(2^{10894}\,11568213483477\cdot10^{18600} > 2^{45}\,11421465^{3100}\). The 3,100-slot run took 135.18 seconds and 638,776 KiB maximum process RSS with six threads on the local Intel. This time includes constructing the graph and checking its first ten unrestricted root counts. [Recorded transfers](../../../output/research73/growth_return/return26.txt) · [Resource log](../../../output/research73/growth_return/return26-time.txt)

Let \(M_n\) count the ordinary connected meanders on \(2n\) labelled sites, with the two banks ordered, and put \(\mu=\lim_{n\to\infty}M_n^{1/n}\). The transfer follows the pairing of the ends crossing a cut and retains their physical order. Let \(o\) be the unique state immediately after the first site: two paired ends, one on either bank. We impose the width cap \(W=26\) at *every* single step. The same component may close only in the final forced step. After deleting the forced initial and final steps, the remaining words are precisely returns to \(o\). Thus a return in \(m\) two-site slots yields one connected meander of order \(m+1\).

Write \(A\) for the two-step transfer and \(R_m=(A^m)_{oo}\). Concatenating \(k\) root returns gives \(R_m^k\) different returns in \(km\) slots, and hence \(M_{km+1}\ge R_m^k\). Taking roots and limits gives \(\mu\ge R_m^{1/m}\), with denominator **\(m\)**. This argument only uses paths that actually return. States unreachable from the root and states unable to return do not affect the count. The existence of \(\mu\) also follows from the usual supermultiplicative seam splice.

## Two integer comparisons

Both start with one nonzero entry \(x_0(o)=S=2^{45}\). In the first, choose an integer shift \(b_i\in\{3,4\}\) before each two-step transfer and set \(x_i=\lfloor A x_{i-1}/2^{b_i}\rfloor\) entrywise. With \(E_m=\sum b_i\), positivity gives

\[
 2^{E_m}x_m\le S A^m e_o,\qquad
 (p/q)^m<R_m\quad\text{if}\quad S p^m<2^{E_m} x_m(o)q^m .
\]

The C++ program computes the two constituent single steps separately. All physical split positions are present, with left and right distinct, and every permitted move contributes once. To obtain the constants in the check file the Python verifier compares integer powers directly. The mantissa threshold bounds \(x_i\) by \(2^{47}\) at each step. A single-step row has at most four entries, so the unscaled two-step sum is at most \(2^{51}\).

For an independent check of the normalization we also iterate \(x_i=\min(2^{50}-1,\lfloor 5 A x_{i-1}/57\rfloor)\) entrywise. Saturating and rounding only lower the result. At \(m=1859\) the root first reaches 35,323,046,266,550, strictly above \(S=35,184,372,088,832\), so \(R_{1859}>(57/5)^{1859}\). No coordinate saturated in this calculation. The two experiments use the same physical transfer. The spectral-vector comparison in the companion branch generates its graph separately and imposes a different final inequality. [Fixed-ratio run](../../../output/research73/growth_return/fixed26.txt) · [Verification](../../../output/research73/growth_return/verified.json)

The generator enumerates matchings by the partner of the first point, its interior matching and its remaining matching. An opening inserts two new partners at the slit. A closing joins the two old partners across the slit, unless they were already paired; the latter case would close a component prematurely. Moving the slit leaves the pairing unchanged. Bracket flags from the full array of partners are lookup keys. No graph component is selected and reflected positions retain separate indices. The small independent check lists full pairs of upper and lower matchings, follows alternating edges to test connectedness, and measures their maximum width at each physical cut.

## Reproduce

Run from the repository root on a C++17 compiler with OpenMP. The transfer is integer and the number of workers does not enter any inequality. The archive contains the log files read by the Python verifier.

```sh
g++ -O3 -std=c++17 -fopenmp research/review73/growth_return/return.cpp -o /tmp/root-return
D=output/research73/growth_return
for w in 2 4 6 8 10 12; do /tmp/root-return "$w" 2 binary 1 1 > "$D/small$w.txt"; done
/tmp/root-return 26 3100 binary 6 100 > "$D/return26.txt"
/tmp/root-return 26 3100 11.4 6 100 > "$D/fixed26.txt"
python research/review73/growth_return/verify.py --write
python research/review73/growth_return/verify.py
```

The archived small logs may use platform-specific line endings; their numerical records use the same format. The SHA-256 values in the check file identify the exact bytes used in that run.
