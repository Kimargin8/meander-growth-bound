# Changes from version 1.0.0

The original width-32 strip and the independent one-root calculation at width 26 are retained. The revised paper adds the following results for ordinary connected closed meanders.

* Admitting excursions to width 34 yields the exact lower bound 12401248195833/1058038773811 for the growth constant. The old strip bound was 4581442652581/392062763189. Two independently indexed programs check the new witness on 386,058,577 odd coordinates.
* For every positive order n, the ratio of the exact count to the new lower bound raised to power n is at least 8.716941 times 10 to the power −8. This combines the new positive vector with the original finite return certificate.
* The deterministic prescribed-relative-error calculation has an upper bound of O*(2^(14n/25) ε^(−3/2)) on bit operations and bits of memory. The previous exponent in n was 4/7. The finite width required by a particular guarantee depends on its prefactor; the code compares the available guarantees.

The supplementary archive includes the integer data, C++ and Python sources and reproduction instructions. The companion paper is included in its original version.

