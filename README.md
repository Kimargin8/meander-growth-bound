# A certified lower bound for the growth constant of plane meanders

**Iaroslav Argunov · November 2026**

Let $M_n$ count the ordered pairs of noncrossing matchings on $2n$ sites whose union is one loop. The paper proves that their exponential growth constant $\mu=\lim_{n\to\infty}M_n^{1/n}$ satisfies

$$
\mu \geq \frac{4581442652581}{392062763189} > 11.6854827408 .
$$

The lower bound comes from a finite transfer with at most 32 open ends. A second construction counts returns to one state in a smaller transfer and proves $\mu>11.421465$. The article also gives a deterministic relative approximation to $M_n$ at prescribed accuracy.

| Manuscript | Files |
|---|---|
| A certified lower bound for the growth constant of plane meanders | [Paper](growth/output/research73/paper/growth_constant.pdf) · [LaTeX source](growth/research/review73/star/growth_constant.tex) |
| Shortest and near-shortest halves of closed plane meanders | [Companion paper](companion/shortest_meander_halves.pdf) · [LaTeX source](companion/shortest_meander_halves.tex) |

## Verify the result

Run the following command from this directory with Python 3.10 or later, using the standard library.

```sh
python3 -B verify.py
```

This checks the distributed files and integer certificates, traces small meanders directly, and reproduces the companion calculations. The concluding line reports success. The largest transfer is a separate native calculation. Rebuilding it requires a 64-bit Linux system with GCC and OpenMP. The independent vector check and the programs for recomputing positive-vector bounds are described in the [growth verification instructions](growth/README.md). The [companion guide](companion/README.md) identifies the program for each statement about cuts.

The recorded [exact bound](growth/output/research73/growth/growth-bounds.json), [independent vector check](growth/output/research73/growth/vector32-final.txt) and [return certificate](growth/output/research73/growth_return/verified.json) can also be inspected individually.

## Citation and reuse

Cite Iaroslav Argunov, *A certified lower bound for the growth constant of plane meanders* (2026). For a result specific to shortest halves, cite the companion paper by its own title. The original code is available under the MIT License. The original manuscripts and certificate data are available under CC BY 4.0. The external reference table retains its OEIS attribution and license. [Citation metadata](CITATION.cff) · [License scope and third-party credit](NOTICE.md)

Archived preprint and supplement: [10.5281/zenodo.24448957](https://doi.org/10.5281/zenodo.24448957).
