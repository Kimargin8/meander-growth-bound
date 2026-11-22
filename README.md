# A certified lower bound for the growth constant of plane meanders

**Iaroslav Argunov · Revised manuscript and verification materials · Version 1.1.0**

The revised paper proves that the growth constant of connected closed meanders is at least $12401248195833/1058038773811>11.7209770594$. The certificate sums selected walks outside a connectivity strip and uses two distinct additional weights on its boundary shifts. For prescribed relative error $\epsilon$, the paper obtains a deterministic algorithm using $O^*(2^{14n/25}\epsilon^{-3/2})$ bit operations and bits of memory.

[Revised paper](paper/growth-bound.pdf) · [LaTeX source](paper/growth-bound.tex) · [Companion paper](companion/shortest_meander_halves.pdf) · [Zenodo v1.1.0](https://doi.org/10.5281/zenodo.24449944) · [First release](https://doi.org/10.5281/zenodo.24448957)

## Reproduce the revision

From the root of this release, using Python 3.10 or later, run

~~~text
python -B verify.py
~~~

The command checks the release manifest and integer comparisons, recomputes the boundary weights in two ways, and traces small diagrams. To reproduce the finite-width choices or run a small confined count, use

~~~text
python -B extension/research/review75/paper/precision.py 1000 1/100
python -B extension/research/review75/paper/precision.py 5 1/100 --count
~~~

The optional large computation uses 64-bit little-endian Linux with GCC and OpenMP. Choose a writable directory on a disk with at least 10 GB free. From the **extension** directory use the following commands with a fresh working directory. The first reconstructs the two earlier vectors, and the next five construct and independently verify the final vector.

~~~sh
python3 -B research/review74/boundary/replay.py --full --work /tmp/meander-v11 --threads 12
python3 -B research/review75/growth/build.py /tmp/meander-v11/gap
OMP_NUM_THREADS=12 /tmp/meander-v11/gap/compute /tmp/meander-v11/short-exits-32.bin 80 /tmp/meander-v11/final.bin
OMP_NUM_THREADS=12 /tmp/meander-v11/gap/audit /tmp/meander-v11/final.bin
g++ -O3 -std=c++17 research/review75/growth/compare.cpp -o /tmp/meander-v11/gap/compare
/tmp/meander-v11/gap/compare /tmp/meander-v11/baseline-32.bin /tmp/meander-v11/final.bin
sha256sum /tmp/meander-v11/final.bin
~~~

The final vector is 3,088,468,656 bytes. Its SHA-256 digest is 23fdd8c621fedbb2ef40556a1bf484758149ed28afebee0aed0a0675123ef538. The [verification receipt](extension/output/research75/growth/verified.json) gives the reduced fraction, the older input digests, the finite lower bound and the sufficient widths. The [derivation](extension/research/review75/growth/reading.md) explains the coordinate comparisons. Python with NumPy and SciPy is needed only for the separate small-width spectral experiments supplied with the earlier excursion certificates.

The [original strip study](growth/README.md) and [companion](companion/README.md) retain their own checks. Run python -B verify_baseline.py from the release root to run the original suite, including its companion calculations. The exact enumeration benchmark on the Intel workstation is a separate project experiment.

See [version changes](CHANGES.md), [citation metadata](CITATION.cff) and [licenses and attribution](NOTICE.md).

