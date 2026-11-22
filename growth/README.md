# A lower bound for the growth constant of closed meanders

**Iaroslav Argunov.** Computational supplement for independent reading.

The ordered upper and lower noncrossing matchings occupy \(2n\) labelled sites, and their union has one component. The width-32 certificates prove **\(\mu>11.685482\)** for their unrestricted exponential growth per order \(n\). On the strip, the best certified lower ratio is \(18325770610324/1568251052756\) and the best upper ratio is \(15264493432827/1305926315603\). The two extrema come from different positive integer vectors. Confining the total number of exposed ends gives the unrestricted lower bound. A separate construction counts returns to one physical state and independently proves \(\mu>11.421465\).

Read the [manuscript](output/research73/paper/growth_constant.pdf) or its [LaTeX source](research/review73/star/growth_constant.tex). The [spectral certificates](output/research73/growth/growth-bounds.json) give exact interval endpoints and state counts. The [primary width-32 transcript](output/research73/growth/i32.txt) records the numerators, denominators and extremal rows of the ordinary comparison at every reported iteration. The [root-return certificate](output/research73/growth_return/verified.json) records the two integer normalization experiments and their direct checks.

To typeset the manuscript separately with a TeX Live installation, run `latexmk -pdf -outdir=build research/review73/star/growth_constant.tex` from the extracted top directory. The bibliography is included in the source.

## Check the supplied evidence

Extract the entire ZIP and run from its top directory with Python 3.10 or newer. The audit uses the standard library.

```sh
python3 verify_package.py
```

The audit checks all manifest hashes and verifies the numerical comparisons with integer arithmetic. It reconstructs small meanders by gluing complete pairs of noncrossing matchings, checks the displayed strip counts, intersects all recorded rational intervals for the same width, and verifies the prescribed-accuracy calculation. The last line prints the lower endpoint for the widest certified strip, rounded down. Large-width row extrema in the transcripts come from the C++ computations. To check those extrema against every row of the actual transfer, regenerate the corresponding run.

## Recompute the transfers

A 64-bit Linux environment with GCC supporting C++17, unsigned 128-bit integers, and OpenMP suffices. The thread count can be changed. Compiled binaries and vectors are temporary.

```sh
python3 research/review73/growth/verify.py --small-native --threads 4
python3 research/review73/growth/verify.py --replay w26 --replay f28 --threads 6
python3 research/review73/growth/verify.py --replay f30 --threads 12
python3 research/review73/growth/verify.py --replay i32 --threads 14
```

The first line compiles three transfers and compares their answers with direct diagrams at seven widths and all orders through eight. The next line regenerates two complete certificates, one with all physical states at width 26 and one with the reflected representation at width 28. The remaining commands reproduce the two larger widths. The archived [clean-build record](output/research73/growth/clean-replay.json) identifies which large runs were repeated from a fresh extraction, the compiler, and the input hashes. The [implementation notes](research/review73/growth/README.md) explain the slit, the forced first and last moves, and the root-accessibility tests.

For a different graph construction and proof, reproduce the [return experiment](research/review73/growth_return/reading.md). It builds all physical partner arrays at width 26, spells out the four legal moves, and derives its inequality by concatenating returns to a single state. It checks 42 width-restricted counts by independent enumeration of complete small diagrams. The 3,100-slot calculation used 624 MiB maximum process RSS on the machine of the original run.

## Choose a strip for a prescribed relative error

The second result gives an explicit sufficient strip width for an order and a rational tolerance \(0<\varepsilon<1\). Counting exactly in that strip gives a deterministic relative approximation to the unrestricted \(M_n\). The bit-operation upper bound is \(O^*(2^{4n/7}\varepsilon^{-3/2})\), where polynomial factors in the input parameters are suppressed as specified in the paper. For a finite example run

```sh
python3 research/review73/growth/accuracy.py --order 1000 --epsilon 1/100
python3 research/review73/growth/accuracy.py --order 6 --epsilon 1/100 --count
```

The first command returns **628**, the least even width meeting the program's explicit sufficient bound in this example. The second command also runs the exact strip transfer and returns 1,828 at order six. The [precision examples](output/research73/growth/accuracy.json) compare `elementary_width`, based on the uniform \(2^{-52}\) factor, with `certified_width`, which also uses \(M_n\ge 2^{-21}\,11.68^n\) at every order. At order 1,000 and one-percent tolerance these widths are 762 and 628. The [avoidance program](research/review73/star/approximation.py), [rational comparisons](output/research73/star/approximation.json), [width-32 edge audit](output/research73/growth/edges32.txt) and [finite prefactor transcript](output/research73/growth/finite32.txt) supply the underlying bounds.

For the full width-32 edge audit, run `python3 research/review73/growth/verify.py --edges32 --threads 14`. To regenerate the finite prefactor, run `python3 research/review73/growth/verify.py --finite32 --threads 14`. The latter computes the integer vector at round 300 with `implicit_bound 32 300 witness.bin`, then invokes `finite_bound witness.bin 40` on the reconstructed data. The run generates its own checkpoint.

The [independent vector verifier](research/review73/growth/verify_vector.cpp) takes the saved odd-slit vector and reconstructs the width-32 componentwise inequalities. The generator saves the last multiple of 100 reached; the round is part of the checkpoint header. To keep the checkpoint and run the independent calculation explicitly:

```sh
g++ -std=c++17 -O3 -fopenmp research/review73/growth/implicit_bound.cpp -o /tmp/strip-vector
g++ -std=c++17 -O3 -fopenmp research/review73/growth/verify_vector.cpp -o /tmp/check-vector
OMP_NUM_THREADS=6 /tmp/strip-vector 32 300 width32.checkpoint
OMP_NUM_THREADS=6 /tmp/check-vector width32.checkpoint --shift
```

The final [independent check](output/research73/growth/vector32-final.txt) reports its round, the dimensions of both slit parities, the root value and the two extremal ratios. The package audit checks both ratios against the corresponding round in the primary transcript. The final two lines record the separate positivity checks for \((A-23I/2)v\) and \((A-58I/5)v\). The audit intersects each available interval with the ordinary one.

The closed-meander values used to check the orders with unrestricted width come from OEIS A005315. The finite-frontier method and the published lower bound used for comparison are due, respectively, to Jensen and to Albert and Paterson. [OEIS A005315](https://oeis.org/A005315) · [Jensen](https://doi.org/10.1088/0305-4470/33/34/301) · [Albert and Paterson](https://doi.org/10.1016/j.jcta.2005.02.006)
