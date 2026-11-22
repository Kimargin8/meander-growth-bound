# A finite lower bound from unfinished outer excursions

Let the two banks and the linear order of sites be distinguished. Start after the compulsory opening site in the one-pair root $o$. The root has odd slit position. The ordinary width-32 transfer and its bank-folded representation are those of [the first preprint](https://doi.org/10.5281/zenodo.24448957). Each odd stored coordinate admits an ordinary path to $o$ in exactly 30 sites, as verified there. A legal shift followed by such a path takes any even coordinate to $o$ in 31 sites. Ordinary two-site returns at $o$ allow padding.

For each directed old boundary shift from $u$ to $u+d$, $d=\pm1$, add a private copy of the path $0,\ldots,34$. The opening from its source enters at $a=u+1$. Keep both directions of every edge in the copy. Add just one exit, from $b=a+d$ to the old destination $j$. Each edge is a single labelled site. Write $P$ for the path adjacency and $A$ for the expanded transfer (source rows). The ends at height 34 retain their pairing until a close. A close at $b$ joins one new end to its old neighbour and restores the original pairing at slit $u+d$. For a complete excursion its first return is fixed by the heights; its closing gap fixes $d$, and its shifts recover the walk in the private path. Hence different expanded complete returns spell different physical words. At the end of an unfinished excursion two guesses for $d$ may still be possible. Both claims remain true if each private path is restricted to an interval containing $a,b$.

## A subeigenvector on the finite expanded graph

Fix $s=\sqrt X>1+\sqrt5$. Suppose $C$ is the stored nonnegative one-site comparison matrix including the ordinary edges and boundary bonuses at most $G_{ab}$, where $G=(sI-P)^{-1}$. The comparison matrix must count both physical shift multiplicities and have row sums at most four. Suppose $C^2v\ge Xv$ on the odd class, $v>0$. Define $f_i=v_i$ on old odd states and $f_i=(Cv)_i/s$ on the old even states (only the even-to-odd block occurs here). On the private copy exiting to $j$ put $f_k=G_{kb}f_j$. There $(sI-P)f=e_b f_j$. On the old states the entry contribution from the opening is $G_{ab}f_j$. Thus $Af\ge sf$ everywhere. These formulas apply to the actual stored coefficients, including dyadic downward comparisons.

Write $V=\max_i v_i$. The infinite two-sided path has diagonal Green function $1/\sqrt{s^2-4}$, since its length-$2h$ returns are counted by $\binom{2h}h$. Embedding the finite path in the line bounds each finite diagonal entry. The matrix $G$ is positive definite for $s>2$, so $G_{kb}\le\sqrt{G_{kk}G_{bb}}\le1/\sqrt{s^2-4}$. Hence

$$
\max_k f_k \ \le\ \frac{\max(1,4/s)}{\sqrt{s^2-4}}\,V
\ \le\ V/2
$$

on all private vertices. The last two sufficient inequalities are $X\ge8$ and $X(X-4)\ge64$, both valid for $s>1+\sqrt5$. The same proof works for a truncated private interval with its own Green function.

## Completion in 32 physical sites

The ordinary reachability test transfers to the physical strip: the root indicator is invariant under exchanging the two banks, so applying 30 ordinary steps to it gives the same count on each physical state as on the corresponding stored coordinate.

Let $m$ be even. Starting at $o$, all endpoints at time $m$ have odd slit position. Give an old vertex terminal weight 1 and a private vertex terminal weight $1/2$, forming a vector $g$. On these odd endpoints, $f/V\le g$. Terminal weights are bank invariant. Consequently the stored value $(A^m g)_o$ equals the corresponding weighted sum of paths in the full physical annotated graph. One way to check the equality is to reflect a state, its chosen exit (if present), and all subsequent sites together. On each half of the slit range the recurrence agrees with the stored recurrence. At the middle the two shift directions keep their labelled multiplicities. This argument uses the invariant terminal vector $g$; no symmetry assumption on a proposed $v$ is needed.

A physical site word ending in the old strip has at most one complete annotation. A word ending outside has at most two, differing only in the guess for its last, unfinished excursion. The sum of its terminal weights is therefore at most one. The number of distinct physical prefixes of length $m$ in this projection is at least $(A^m g)_o\ge s^m v_o/V$.

Every such prefix has a completion to $o$ in exactly 32 further sites in the unrestricted physical language. An old odd endpoint needs 30 plus a two-site return. A private endpoint has odd slit $k$, so $1\le k\le33$. If the two ends adjoining this slit have different partners, close them immediately. The old destination has even slit and takes 31 more sites. If the ends are paired together, choose an adjacent internal slit, shift once and close there. The new tested pair shares one end with the former pair and cannot itself be matched. The destination has odd slit and takes 30 more sites. The possible final close here is unrestricted; an unfinished announced direction need not be honoured. All appended paths are physical and have no premature closed component. Selecting one completion for each distinct prefix gives an injection into ordinary connected meanders: the first $m$ interior sites recover the prefix.

With the two compulsory sites this proves, for every $n\ge17$,

$$
M_n\ \ge\ (v_o/V)X^{\,n-17}.
$$

Thus $M_n\ge 2^{-p}X^n$ at all positive orders whenever $2^p v_o\ge V X^{17}$ and $2^p\ge X^{16}$ (use $M_n\ge1$ up to order 16). The number 32 cannot be reduced for a uniform completion of all odd physical endpoints. From height 34 and slit 1 to the root, if $c,a$ count closes and opens, $c-a=16$. Net upper-opening slit shifts must equal $c-a$ as well, so the total length is at least 32. The height-34 slit-1 position is reachable by opening from a width-32 slit-0 state. For example, reach width 32 by successive openings and shift to slit 0.

If completion is required to respect the private annotation, the original uniform bound is 64 sites. The distance to the internal target is at most 33, plus one exit and at most 31 original sites, a total of at most 65. The whole path has even length and therefore costs at most 64. With restricted copies $|k-b|\le4$ the same proof gives 36. The projection used above yields 32 for both models. These completion statements are distinct from the equality of complete returns.

## Values from the width-32 witness

The round-440 independent scan has $X=12401248195833/1058038773811$, $v_o=245357193322677$ and $V=2667940370280775$. Exact comparisons give $p=64$ and $M_n\ge2^{-64}X^n$ at all orders. With $Q=159/100$ and $r=3898274927/10^9$ from the published height estimate, $(r^2/X)^{25}<Q^{14}$ and $Q^3>4$. Section 5 of the first preprint bounds the count omitted at width $w$ by $2(2n-1)r^{2n}Q^{-w}$. Taking the least sufficient even width, capped at full width, now proves prescribed relative error $\epsilon$ with $O^*(2^{14n/25}\epsilon^{-3/2})$ bit operations and bits of memory.

The exact input values and thresholds can be checked without third-party dependencies:

~~~text
python -B research/review75/proof/check_prefix.py --x 12401248195833/1058038773811 --root 245357193322677 --maximum 2667940370280775
~~~
