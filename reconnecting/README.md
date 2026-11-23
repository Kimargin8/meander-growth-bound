# Reconnecting first returns at width 32

The four files in [certificates](certificates) ending in `-audit.txt` or `-rescan.txt` report two successive readings of the positive vectors. Each reading checks every coordinate. The mixed comparison rounds down and uses the minimum. The unit comparison rounds up and uses the maximum. The growth constant is per pair of crossings, the river has an ordered set of sites and the two banks are distinguished.

| Comparison | Inequality | SHA-256 of the final vector |
| --- | --- | --- |
| Reconnecting returns | Lower, 1461549008357498 / 124174541777295 | 192394a2f012fcc6e9913767c0916ec3cc71b1f977da7a3fb5c1f8d4d339f334 |
| Nested unit returns of every depth | Upper, 3431031596329 / 292144566372 | 0afb1f908e93533bcd7493be3758b4c0aaafb16720d7233a56ee06ca384b90ca |

## Small checks

From the release root run `python -B verify.py` with Python 3.10 or later. The [lightweight checker](certificate.py) compares the lower weights from path minors with a positive integer series through length 83. For the upper comparison it checks two upper majorants and the infinite-line cap. For each of three strip widths, three nesting depths and seven orders it compares two ways of counting complete words. It also checks the extrema reported in the large scans and the rational comparisons in the approximation theorem. This command uses the Python standard library. The multi-gigabyte witness vectors are reproducible from source and are omitted from the archive.

## Reconstruct both width-32 comparisons

Use little-endian 64-bit Linux, Python 3.10+, GCC with OpenMP, approximately 20 GiB available RAM and 25 GB of free disk space. From the release root, with a freshly chosen working directory, run

~~~sh
python3 -B reconnecting/replay.py --work /tmp/meander-v12 --threads 8
~~~

The driver reconstructs rounds 300, 360 and 440 using the earlier strip programs, then round 500 and the distinct mixed and unit witnesses at rounds 560 and 540. It builds and runs the separate physical-coordinate readers and compares the vector digests, coordinate counts and rational extrema. The generator connects two old partners before inserting a new pair. The mixed reader reconstructs the 34-end diagram and performs the splice on its actual partners. The units for memory above cover the original construction of the mixed witness.

If a round-440 vector from version 1.1 is already available, supply `--from-v11 /path/to/round440.bin`. To audit saved vectors directly, provide round 500 followed by both final inputs:

~~~sh
python3 -B reconnecting/replay.py --work /tmp/meander-audit --threads 8 --audit-vectors /path/to/round500.bin /path/to/mixed560.bin /path/to/unit540.bin
~~~

The driver verifies the round-500 SHA-256 (d7bd9f5549906c7bf5a7e1d591cd083e4c33a511d51369802c66cfb8c03802a2). Each vector is 3,088,468,656 bytes including its header. A Linux C++20 compiler is needed for the smaller, marked local-window calculations in the approximation theorem. Run

~~~sh
python3 -B reconnecting/replay.py --work /tmp/meander-window --threads 8 --windows-only
~~~

This regenerates the positive vectors and independently decodes and checks the successor lists for all 716,861 and 2,733,515 states. Source files are under [native](native) and all four certificates use integer inequalities. The upper bound for the unit-return language at width 32 and the separate upper bound for all plane meanders have different scopes. The unrestricted upper bound is proved in [Meandric growth from local connectivity at a cut](https://doi.org/10.5281/zenodo.24457226).
