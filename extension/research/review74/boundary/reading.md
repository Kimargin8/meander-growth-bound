# Eliminating the outer layer of a meander strip

Iaroslav Argunov · Research note and reproduction guide · November 22, 2026

The two banks and the order of the sites are distinguished. Let $T_W$ be the one-site transfer on physical path pairings of at most $W$ ends and a split between the two banks. A double close is legal if its ends belong to different paths. Removing the compulsory first and last sites of a connected meander leaves a return to the one-pair root. A step of $T_W$ contributes one site. Two sites increase the order by one.

## An exact formula for excursions of every duration

Put $N=W+2$. Order the states of $T_N$ as the inner states followed by those at height $N$. The off-diagonal blocks are $B$ for opening and $C$ for closing. A move that stays outside only shifts the split, retaining the pairing. The exterior block $D$ is consequently a direct sum of copies of the path adjacency $P$ on $0,1,\ldots,N$, one per outer pairing. At a per-site parameter $s>2$, eliminating the outer block gives the exact operator

$$
S(s)=T_W+B(sI-D)^{-1}C .
$$

The first-return generating matrix is $S(s)/s$. An excursion with $j$ outer shifts contributes $s^{-j-1}BD^j C$ to $S(s)$ and uses $j+2$ sites. Visits to the inner layer uniquely separate first returns. Because $B$ and $C$ have row sums at most one and $D$ has row sums at most two, discarding durations above $h\ge2$ changes $S(s)$ in maximum row-sum norm by at most $s^{-1}(2/s)^{h-1}/(1-2/s)$.

For $s=p/q$, define $K_0=1,\ K_1=p,\ K_{i+1}=pK_i-q^2K_{i-1}$. With $m=N+1$ and indices $0\le i,j<m$, the common-denominator kernel is

$$
K_m[(sI-P)^{-1}]_{ij}=q^{1+|i-j|}K_{\min(i,j)}K_{m-1-\max(i,j)}.
$$

Substitution in the tridiagonal equations proves the identity, with both exterior boundary values set to zero. For one pairing, assemble the vector of allowed closes, perform forward elimination and backward substitution, and return the result to the inner boundary. Each adjacent paired arc specifies one entrance; deleting it recovers the unique source. This uses $O(W)$ arithmetic operations and $O(W)$ temporary entries per pairing. For rational $s$ the integer bit lengths grow with $W$. For $(W,N)=(32,34)$ the formula still visits 129,644,790 pairings with 35 gaps each. A full pass of this formula at width 32 remains a separate timing experiment.

At $s=7/2$, streamed_path.py checks the exact residual for every outer state at $W=2,4,\ldots,14$. A separate row-based computation inserts the adjacent pair, follows the two former partners and evaluates the cofactor formula. It confirms every returned numerator, including all 6,435 at $W=14$. The program literal.py additionally enumerates complete pairs of arch matchings and traces each cycle. It checks 175 integer counts at prescribed excursion cutoffs.

If the exit returns to one of the two neighbouring gaps, all durations can also be summed without visiting the outer pairings. For $m=W+3$, $D_0=1$, $D_1=s$ and $D_k=sD_{k-1}-D_{k-2}$, the added weight between outside gaps $i,i+1$ is $D_iD_{m-i-2}/D_m$. It depends only on the gap. This is the adjacent entry of the same resolvent. At $W=32$ and $s^2=12$ the 32 weights range from 0.11122565 to 0.11237244, rounded outward. The script adjacent_weights.py independently checks 272 edges by solving path systems at $s=4$. These all-duration weights are a candidate for a subsequent width-32 certificate. The 36-site return argument below uses the finite three-length model.

## Why the restricted language proves a larger lower bound

The accompanying [mathematical addendum](../../../output/research74/excursions/boundary_excursions.pdf) proves that an outer excursion of length 3, 5 or 7 whose shifts end one position from their start can realize the same inner destination as a single legal top shift. At least 1, 3 and 9 such paths, respectively, are available for each direction. The three lengths are distinct. Within a length, different slit walks yield different words.

Let $U$ contain these original top shifts, counted with multiplicities. For a target $L\le12$ per order, the exterior terms contribute at least $1/12+3/12^2+9/12^3=7/64$ times $U$ to $S(\sqrt L)$. Set $H=T_{32}+(7/64)U$. Round the additional nonnegative term down on each sweep. Two sweeps on a positive vector $v$ give an integer lower approximation to $H^2v$. Let $L$ be the smallest observed ratio to $v$. Every stored odd state reaches the root in exactly 30 ordinary sites. The addendum uses these paths to bound root returns of $H$ from below and prove that the root-return series of the admitted physical words diverges at $1/\sqrt L$. Inclusion in the connected language gives $\mu\ge L$.

The independent program reads the proposed integer vector, reconstructs each move on physical ends and verifies the minimum over all 386,058,577 odd coordinates. Both implementations obtain $5660040998411745/483083866680971>11.7164769697$. The reported maximum ratio of the rounded sweep is diagnostic.

## Relative approximation

Write $X=5660040998411745/483083866680971$. Represent each admitted excursion by a path with one edge per site. Every odd vertex in this larger graph reaches the root in exactly 36 sites, including the vertices inside a represented excursion. Section 4 of the addendum proves $M_n\ge 2^{-71}X^n$ at every order. The earlier height estimate then suffices for deterministic relative error $\epsilon$ in $O^*(2^{9n/16}\epsilon^{-3/2})$ bit operations and bits of memory. The earlier exponent was $4/7$. The stronger asymptotic estimate has a smaller finite prefactor. For example, at one percent relative error and order 1000, combining the two proofs selects the earlier width of 626. At order 20000 it selects the new width of 11358. These are sufficient widths computed by integer comparisons.

## Reproduction

Use Python 3.10 or later. Run the commands from the root of the project or extracted packet. The first verifies the recorded arithmetic and exact outer-path calculation. The native tests require 64-bit little-endian Linux, GCC with OpenMP, and a writable working directory. Full regeneration creates two checkpoints of about 3.09 GB each. The additional 60 rounds used 8.54 GiB peak resident memory. The original 300-round construction is included in the full replay.

~~~text
python -B research/review74/boundary/replay.py
python -B research/review74/boundary/relative.py
python -B research/review74/boundary/relative.py --order 1000 --epsilon 1/100
python -B research/review74/boundary/relative.py --order 5 --epsilon 1/100 --count
python -B research/review74/boundary/replay.py --small --work /tmp/strip-small
python -B research/review74/boundary/replay.py --full --work /tmp/strip-full --threads 12
~~~

The independent complete-diagram check and the small-width interval checks also use NumPy and SciPy.

~~~text
python -B research/review74/boundary/literal.py
python -B research/review74/excursions/scan.py --verify
python -B research/review74/excursions/shift_witness.py --verify
~~~

The base transfer and starting vector construction are documented in the [previous preprint and code release](https://doi.org/10.5281/zenodo.24448957).
