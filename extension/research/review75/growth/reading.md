# Adjacent-gap excursions at width 32

The order of the sites and the two banks are distinguished. An interior site either opens, splices or shifts the slit. Write \(T\) for the ordinary one-site transfer with at most 32 exposed ends, and \(o\) for its one-pair root. A closed connected meander has a first opening, a path returning to \(o\), and a last closing. We use the bank-folded storage and the physical transition interpretation in the [original preprint](https://doi.org/10.5281/zenodo.24448957).

## Exact weights

When the width is 32, allow an opening to 34 ends, an arbitrary sequence of outside shifts, and then a splice returning to the neighbouring gap of a formerly legal top shift. The matching remains constant while outside. For a left or right shift from \(u\), the initial outside gap is \(i=u+1\) and the closing gap is \(j=i\pm1\). At this gap the new end splices to the old neighbouring end and the other new end replaces it. This gives exactly the previous inner destination. Height records the departure and first return. In a complete physical word these positions determine its decomposition, direction and outside walk.

Let \(P_m\) be the path on \(m=35\) vertices labelled \(0,\ldots,34\). At per-site parameter \(s>2\), the physical excursion contributes the following extra coefficient to the effective one-site transfer in direction \(i\to j\).

\[
 G_{ij}(s)=[(sI-P_m)^{-1}]_{ij}
          =\sum_{k\ge0}(P_m^k)_{ij}s^{-k-1}.
\]

The opening and closing explain the exponent \(k+1\) after the common one-site factor has been removed. Put \(p=1173/100\), and choose the integer \(q_{ij}=\lfloor4096G_{ij}(\sqrt p)\rfloor\). Among the 32 internal undirected edges these integers equal 469 at the two ends and 474 on the remaining 30. Both orientations agree. For \(4<x\le p\) the same lower coefficients are valid at \(s=\sqrt x\).

Two exact derivations produce the table. The first uses continuants, with \(D_0=1,D_1=s,D_k=sD_{k-1}-D_{k-2}\) and \(G_{ij}=D_iD_{m-1-j}/D_m\) for \(i\le j\). The second counts walks one shift at a time and bounds the unexamined tail using at most two choices per shift and the geometric ratio \(4/p\). All 272 internal edges over even widths from 2 to 32 are checked both ways. The weights are invariant under reflection. Root-return coefficients in physical and stored addresses agree by applying the first-return series to the invariant root indicator, as in the preceding addendum.

A smaller language has the same dyadic lower coefficients. Restrict the outside gap to distance at most **four** from its prescribed closing gap. For each directed edge, recompute \(G\) on the intersection of the path with this radius-four interval. The integer floors for both directions of every width-32 edge remain exactly 469 or 474. This fact permits a finite terminal estimate for arbitrarily long outside walks.

## Integer certificate

Add \(q_{ij}/4096\) times each old top shift to the ordinary one-site matrix, obtaining \(C\). Apply downward rounding once to the sum of the extra terms in each row. The resulting monotone integer rule \(F\) satisfies \(F(F(v))\le C^2v\). From the previous round-360 witness, 80 power iterations produce the round-440 positive odd vector \(v\). The generator and physical reconstruction audit compare every odd component after two sweeps. The audit reconstructs arcs directly, computes reflection and splice afresh, stores odd and even arrays separately, and uses the independent walk table. Their common conclusion is

\[
 C^2v \ge Xv,\qquad X=12401248195833/1058038773811 > 11.7209770594 .
\]

Both implementations compare ratios by integer cross-products. With normalized inputs at most \(2^{52}\), row sums below or equal to four give at most \(2^{54}\) at the intermediate sweep. The sum of the products with the two integer weights is less than \(2^{64}\) because each weight is strictly below 512. Final sums are at most \(2^{56}\); ratio comparisons fit unsigned 128-bit arithmetic. The independent sweep uses 128-bit arithmetic for the extra weighted sum as well.

At \(s=\sqrt X\) the physical first-return matrix including excursions dominates \(C/s\) entrywise. Original odd states have ordinary 30-site paths to \(o\). The positive witness hence forces divergence of the admitted root-return series at \(1/s\), by the same root expansion as in the previous addendum. This proves \(\lim M_n^{1/n}\ge X\).

## Uniform finite counts

Represent an excursion by first choosing one of its permissible closing gaps. Use a private outside path restricted to radius four of that closing gap and allow final splicing only there. Partial prefixes of the two choices can coincide physically. Every complete word specifies the actual closing gap at its return, so complete rooted walks retain exactly the physical multiplicities. Folded first-return coefficients agree with the physical ones by reflection. Write \(A\) for the one-site adjacency of this finite stored representation. Give old odd vertices the weights \(f=v\), old even vertices \(f=Cv/s\), and an outside vertex with current gap \(i\), closing gap \(j\), and destination \(t\) the weight \(G_{ij}(s)f_t\) on its restricted interval. The resolvent equation gives \(Af=sf\) at outside vertices. At original vertices \(Af\ge Cf\ge sf\). Every odd vertex reaches the root in 36 sites: at most four shifts to the chosen gap, one closing, and at most 31 ordinary sites from the destination, with two-site root returns for padding.

For a finite constant, reuse the original 80-site return certificate for the old round-300 witness \(v^0\). It gives \(T^{80}e_o\ge 2^{76}(308588629582/208548141853)v^0\) on old odd vertices. Comparing both normalized vector files coordinate by coordinate supplies \(g=\min_i v^0_i/v_i\). Put \(\gamma=2^{76}(308588629582/208548141853)g\). Thus \(A^{80}e_o\ge\gamma f\) at all old odd vertices.

Write \(w=474/4096\). On original vertices \(Tf\ge sf/(1+w)\), since \(C\le(1+w)T\). From an original odd vertex allow six ordinary steps before using the 80-site certificate. From an odd auxiliary vertex at distance \(k-1\) from its target, \(1\le k\le5\), take its shortest path to the target, close, and take \(6-k\) ordinary steps. This ends at an old odd vertex. Consequently \(A^{86}e_o\ge\gamma\beta s^6 f\) on the odd class, where \(\beta\) is the minimum of the original factor \((1+w)^{-6}\) and the auxiliary factors

\[
  \frac{1}{(1+w)^{6-k}s^k G_{ij}(s)} .
\]

The scaled Green entry \(s^kG_{ij}(s)\) is nonincreasing in \(s^2\), by the walk series whose first walk has length \(k-1\). Evaluation at the smaller rational \(s^2=117/10\) gives a uniform certified lower \(\beta\). The checker covers every radius-four interval, direction and current gap. It follows for \(n\ge44\) that

\[
 M_n\ge (A^{2n-2})_{oo} \ge \gamma\,\beta\,v_o X^{n-41}.
\]

For \(1\le n\le43\) use the published bound \(M_n>4.71986\,10^{-7}L_0^n\) with \(L_0=4581442652581/392062763189\). The verifier takes the lesser of the two resulting constants and rounds it down. Combined with the previous height estimate, \((r^2/X)^{25}<Q^{14}\) at \(r=3898274927/10^9\), \(Q=159/100\). Deterministic approximation from below with relative error \(\epsilon\) therefore costs \(O^*(2^{14n/25}\epsilon^{-3/2})\) bit operations and bits of memory. The previous rational exponents were \(4/7\) in the public paper and \(9/16\) at stage 74. The exact count at order 39 remains open in this experiment.

The coordinate check gives \(g=158672161814/348541137730\), \(v_o=245357193322677\). The terminal factor is \(\beta=300444138341140529152/640297620801287196875\). At every order the resulting explicit inequality is

\[
                  M_n \ge 8.716941\cdot10^{-8}\,X^n .
\]

For a relative error of one percent at \(n=1000\), the published certificate requires width 626, and the new one requires 624. At \(n=20000\) the widths are 11398 and 11270. At \(n=39\) they are 76 and 78; the earlier certificate gives the smaller sufficient width there. These examples evaluate sufficient widths by integer comparisons.

## Reproduction

From the project root run the commands below under 64-bit little-endian Linux with Python 3.10+ and GCC. Use distinct writable input and output checkpoint paths. To regenerate the older vector use the preceding verification package and the full-replay instructions in the stage-74 boundary guide. The three checkpoints are each about 3.09 GB. The stage-75 generator used 8.54 GiB of resident memory and 4 min 2.34 s on 12 threads. The independent audit took 40.10 s including reading the witness from the host filesystem. All three digests are included in the verification receipt.

~~~text
python3 -B research/review75/growth/build.py /tmp/growth75
OMP_NUM_THREADS=12 /tmp/growth75/compute data/interim/review74-boundary/width32-round360.bin 80 data/interim/review75-growth/width32-round440.bin
OMP_NUM_THREADS=12 /tmp/growth75/audit data/interim/review75-growth/width32-round440.bin
g++ -O3 -std=c++17 research/review75/growth/compare.cpp -o /tmp/growth75/compare
/tmp/growth75/compare output/research73/growth/width32.local.checkpoint data/interim/review75-growth/width32-round440.bin
python3 -B research/review75/growth/verify.py
~~~
