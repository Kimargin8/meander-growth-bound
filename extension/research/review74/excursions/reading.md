# First returns across the boundary of a meander strip

Let the two banks be labelled and number all exposed ends in their boundary order. A state consists of its noncrossing path pairing and the upper slit position. Every retained physical state at widths from 2 up to the even integer $W$ is included. In this file $\lambda$ measures growth per **site**. The growth constant per order is $x=\lambda^2$.

## Finite experiment

A first return to widths at most $W$ has either length one or the following form. A double opening at width $W$ raises the width to $W+2$. While outside, only single-bank shifts are permitted. A double close returns to width $W$ if its two ends belong to different paths. Inside and outside an interval the usual prohibition against premature closed components still applies. The position of each successive first return is read directly from the height sequence, so every admitted labelled word decomposes uniquely. The experiment excludes any visit above $W+2$.

Let $T$ be the single-site transfer on the original strip and $E_d$ the matrix of first returns of length $d>1$. There is no length-two return. With a cap $h$ on excursion length put

$$
B_h(z)=zT+\sum_{3\le d\le h} z^d E_d.
$$

The root-return generating function is the root entry of $(I-B_h(z))^{-1}$. For $z>0$ the ordinary transitions remain irreducible and $B_h(z)$ is strictly increasing in spectral radius. Consequently the per-site growth rate is the positive $\lambda$ at which $\rho(B_h(1/\lambda))=1$. This also follows by replacing a weighted edge of length $d$ with a path of $d$ unit edges. In either graph the root returns have even length, and the ordinary two-site return permits concatenations at all subsequent orders. This justifies the square converting per-site bounds to per-order bounds.

The table shows representative decimals. Each endpoint for $\lambda$ used in the certificates has denominator $10^8$ and each bracket has width $3\cdot10^{-8}$ per site. The certified percentage intervals, given as exact fractions in the JSON files, quantify the share of the improvement from the original strip to the full next strip.

| $W$ | Original strip | Returns up to 6 | Full $W+2$ | Share recovered |
|---:|---:|---:|---:|---:|
| 8 | 8.17076321 | 8.93237464 | 9.09190350 | 82.6814% |
| 10 | 9.09190350 | 9.63051305 | 9.74052292 | 83.0394% |
| 12 | 9.74052292 | 10.13237201 | 10.21145571 | 83.2070% |
| 14 | 10.21145571 | 10.50418364 | 10.56295055 | 83.2809% |

For width 14 the fixed-$z$ matrix with length cap six has 45,244 nonzero entries on 8,787 states, compared with 91,070 on 33,097 states for the full width-16 matrix. CSR buffers and two floating-point vectors total 718,672 and 1,754,784 bytes respectively. Creating the four separate coefficient matrices and ordinary matrix took 0.181 s, compared with 0.111 s to build the full width-16 transfer. In the fixed-parameter single-sweep microbenchmark the median times were 38.8 and 110.7 microseconds. These timings measure local sparse applications. The complete process for eigenvector search and exact certification also allocates workspace, including additional blocks when linearizing the recurrence. The small-width measurements do not establish a speedup for computing an exact unrestricted $M_n$.

Single-process peak memory gives a different comparison from buffer sizes. On width 14 the full certificate search and exact test used 66.74 MB for the base, 81.00 MB for the six-step cap, and 79.02 MB for the full next strip. Construction and certification together took 0.080, 0.763 and 0.261 seconds respectively in three fresh processes. A separate pilot evaluated the short-return polynomial at successive fixed parameters without the delayed companion matrix. It gave the same bracket, 0.543 seconds and 83.03 MB peak. These are single-run measurements on the local Python installation. The spectrally smaller short-return model saves coordinates during application; the full next strip was faster to certify on these inputs. The measurements are in `peak-memory.json` and `fixed-parameter-pilot.json`. The separate-process experiment can be repeated with `python -B research/review74/excursions/profile_memory.py`.

## A sublanguage using only existing addresses

Write $S$ for the ordinary single-bank shifts starting at width exactly $W$. An opening followed by one left or right shift and a close reaches the same labelled address as the corresponding single shift. In fact, after an opening insert new adjacent ends $L,R$. A shift by one before a close splices one of the new ends with the neighbouring old end. The other new end takes that old end's place and partner. The external order therefore restores the original pairing, with precisely the displaced slit. This proves $E_3=S$, with both directions counted individually.

More generally a walk outside that ends at the same adjacent gap, whatever its length, has that same final splice. Each legal shift uses adjacent *internal* vertices of the path on positions $0,\ldots,W+2$. Choose once and for all a five-vertex subpath in which the given edge is either $1\!-\!2$ or $2\!-\!3$ in local indexing. Such a choice exists even at $W=2$. The number of walks of length $2k+1$ along this subpath between the relevant vertices is $3^k$. On the odd vertices $1,3$, two steps have matrix with rows $(2,1),(1,2)$, and the adjacency column for the centre is $(1,1)$. Reversal covers the other orientation. The endpoint definition prevents double counting of the two shift directions. Choosing one subpath for each edge prevents counting the same walk twice through overlapping windows.

For $z=1/\lambda$, these words add $z^3 S/(1-3z^2)$ to $zT$. After factoring out $z$, the improvement in weight on each affected shift is $1/(x-3)$. The first three terms alone give $1/x+3/x^2+9/x^3$, which is at least $7/64$ for $3<x\le 12$. Thus in a test at a target $x\le 12$, $T+(7/64)S$ can serve as a simple lower comparison on the original addresses. The spectral radius squared of this constant-weight transfer is not asserted to equal that of the five-vertex sublanguage.

| $W$ | Gain share for fixed $7/64$ | Gain share from three terms | Gain share from full five-vertex subpaths |
|---:|---:|---:|---:|
| 8 | 43.7993% | 70.7881% | 73.8661% |
| 10 | 47.6487% | 68.6120% | 71.1282% |
| 12 | 50.4200% | 67.0478% | 69.2348% |
| 14 | 52.4769% | 65.8774% | 67.8538% |

## Reproduce the certificates

From the repository root use the following commands with the existing Python environment containing NumPy and SciPy.

```text
python -B research/review74/excursions/scan.py --verify
python -B research/review74/excursions/shift_witness.py --verify
```

Removing `--verify` builds new positive integer witnesses and trials. The ARPACK eigenproblem proposes the vectors. Each submitted vector is independently tested at both rational ends using integer inequalities on every coordinate. At denominator $D$, site numerator $a$, and cap $h$, the scaled residual is $\sum_d D^d a^{h-d}E_d v-a^h v$, with $E_1=T$. The lower test finds a strictly positive minimum and the upper test a strictly negative maximum. Matrix sums against the integer vector are bounded explicitly within signed 64-bit arithmetic before use. Powers, denominators, residuals and their comparisons use arbitrary-precision integers. The corresponding checks for the rational five-vertex generating function clear $a^2-3D^2>0$.

All four widths retain separate certificates for caps 3, 4, 5, and 6, the base and the full next strip. For interval calculation let $b,m,f$ be their per-order growth rates for the base, a model and the full strip. For $b<m<f$ the gain $(m-b)/(f-b)$ increases with $m$ and decreases with either $b$ or $f$. The computation uses the corresponding one-sided endpoints. A separate labelled-end audit compares the full actual $E_3$ to $S$, follows the splice independently, and counts all one-, three- and five-shift walks at every allowed slit for widths $2,4,\ldots,14$.

The `output/research74/excursions/` directory contains the integer vectors, their file hashes, both scaled slacks, exact interval endpoints, matrix sizes and machine-dependent timings. Rerunning the proof checks uses only the stored vectors and the combinatorial transfer. Timings and the approximate eigenvalues do not enter the proof.
