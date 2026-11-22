"""Finite first-return polynomials for one extra layer of the meander strip."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import platform
from statistics import median
import sys
from time import perf_counter

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, eigs

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/research74/excursions"
SOURCE = Path(__file__).resolve()
DEN = 100_000_000


@lru_cache(None)
def matchings(pairs):
    if not pairs:
        return ((),)
    result = []
    for inside in range(pairs):
        end = 1 + 2 * inside
        for left in matchings(inside):
            for right in matchings(pairs - inside - 1):
                result.append((end,) + tuple(i + 1 for i in left) + (0,)
                              + tuple(i + end + 1 for i in right))
    return tuple(result)


def open_pair(pairing, gap):
    def moved(i):
        return i if i < gap else i + 2

    return (tuple(moved(i) for i in pairing[:gap]) + (gap + 1, gap)
            + tuple(moved(i) for i in pairing[gap:]))


def splice(pairing, gap):
    if not 0 < gap < len(pairing) or pairing[gap - 1] == gap:
        return None
    left, right = pairing[gap - 1], pairing[gap]
    merged = list(pairing)
    merged[left], merged[right] = right, left
    del merged[gap - 1:gap + 1]
    return tuple(i if i < gap - 1 else i - 2 for i in merged)


def matrices(width, max_length):
    states = [(pairing, g) for k in range(1, width // 2 + 1)
              for pairing in matchings(k) for g in range(2 * k + 1)]
    address = {state: i for i, state in enumerate(states)}
    edges = {d: ([], [], []) for d in [1] + list(range(3, max_length + 1))}

    def append(d, src, destinations):
        rows, cols, counts = edges[d]
        for dst, count in destinations.items():
            rows.append(src)
            cols.append(address[dst])
            counts.append(count)

    for src, (pairing, gap) in enumerate(states):
        ordinary = Counter()
        if gap:
            ordinary[pairing, gap - 1] += 1
        if gap < len(pairing):
            ordinary[pairing, gap + 1] += 1
        closed = splice(pairing, gap)
        if closed is not None:
            ordinary[closed, gap - 1] += 1
        if len(pairing) < width:
            ordinary[open_pair(pairing, gap), gap + 1] += 1
        append(1, src, ordinary)
        if len(pairing) != width or max_length < 3:
            continue
        outside = open_pair(pairing, gap)
        shifts = Counter({gap + 1: 1})
        # Exit opens an adjacent pair. Immediate splicing would close it into a loop.
        assert splice(outside, gap + 1) is None
        for d in range(3, max_length + 1):
            moved = Counter()
            for g, c in shifts.items():
                if g:
                    moved[g - 1] += c
                if g < len(outside):
                    moved[g + 1] += c
            shifts = moved
            returns = Counter()
            for g, c in shifts.items():
                closed = splice(outside, g)
                if closed is not None:
                    returns[closed, g - 1] += c
            append(d, src, returns)

    size = len(states)
    return {d: coo_matrix((np.asarray(counts, dtype=np.int32), (rows, cols)),
                          shape=(size, size)).tocsr()
            for d, (rows, cols, counts) in edges.items()}


def approximate(parts):
    """Linearize the delayed recurrence; only the first block supplies a witness."""
    n = parts[1].shape[0]
    degree = max(parts)
    if degree == 1:
        matrix = parts[1].astype(float)
    else:
        def product(x):
            y = np.empty_like(x)
            y[:n] = parts[1] @ x[:n]
            for d, part in parts.items():
                if d > 1:
                    y[:n] += part @ x[(d - 1)*n:d*n]
            y[n:] = x[:-n]
            return y

        matrix = LinearOperator((n * degree, n * degree), matvec=product, dtype=float)
    vals, vecs = eigs(matrix, k=1, which="LR", v0=np.ones(matrix.shape[0]),
                      ncv=24, tol=2e-13, maxiter=5000)
    r = vals[0].real
    assert abs(vals[0].imag) < 1e-11
    v = vecs[:n, 0].real.copy()
    if v[0] < 0:
        v *= -1
    assert v.min() > 0
    witness = np.rint(v / v.max() * (1 << 48)).astype(np.int64)
    assert min(witness) > 0
    return r, witness


def integer_products(parts, witness):
    assert witness.dtype == np.int64 and np.min(witness) > 0
    for part in parts.values():
        row_mass = int(np.max(np.asarray(part.sum(axis=1))))
        assert row_mass * int(np.max(witness)) < 2**63
    return {d: mat @ witness for d, mat in parts.items()}


def margins(products, witness, a):
    """Integer residual of sum_d (DEN/a)^d A_d v >= v."""
    degree = max(products)
    coefficients = [DEN**d * a**(degree - d) for d in products]
    denominator = a**degree
    arrays = list(products.values())
    least = greatest = None
    for i, v in enumerate(witness):
        residual = sum(c * int(arr[i]) for c, arr in zip(coefficients, arrays)) - denominator * int(v)
        least = residual if least is None else min(least, residual)
        greatest = residual if greatest is None else max(greatest, residual)
    return least, greatest


def check_enclosure(parts, witness, pair):
    lo, hi = pair
    products = integer_products(parts, witness)
    lower_margin = margins(products, witness, lo)[0]
    upper_margin = -margins(products, witness, hi)[1]
    assert min(lower_margin, upper_margin) > 0
    return [str(lower_margin), str(upper_margin)]


def certify(parts):
    r, witness = approximate(parts)
    center = math.floor(r * DEN)
    for pad in (1, 3, 10, 30, 100):
        pair = [center - pad, center + pad + 1]
        try:
            margin = check_enclosure(parts, witness, pair)
        except AssertionError:
            continue
        lower, upper = [Fraction(t, DEN)**2 for t in pair]
        return witness, {"site_numerators": pair, "site_denominator": DEN,
                         "strict_scaled_slack": margin,
                         "order_interval": [str(lower), str(upper)],
                         "order_approx": round(r*r, 10)}
    raise ArithmeticError("could not isolate a positive integer certificate")


def gain_fraction(base, model, full):
    b0, b1 = map(Fraction, base["order_interval"])
    m0, m1 = map(Fraction, model["order_interval"])
    f0, f1 = map(Fraction, full["order_interval"])
    assert b1 < m0 <= m1 < f0
    return [(m0-b1)/(f1-b1), (m1-b0)/(f0-b0)]


def storage(matrix):
    return matrix.indptr.nbytes + matrix.indices.nbytes + matrix.data.nbytes


def sparse_bench(plain, partial, full, root):
    effective = partial[1].astype(float) / root
    for d, mat in partial.items():
        if d != 1:
            effective += mat / root**d
    matrices = {"base": plain.astype(float) / root,
                "short6": effective,
                "full": full.astype(float) / root}
    result = {}
    for name, mat in matrices.items():
        result[name] = {"states": mat.shape[0], "nnz": mat.nnz,
                        "csr_bytes_float64_int32": storage(mat),
                        "two_vector_bytes": 2 * 8 * mat.shape[0],
                        "matvec_median_us": None}
    times = {key: [] for key in matrices}
    inputs = {key: np.linspace(1., 2., mat.shape[0]) for key, mat in matrices.items()}
    names = list(matrices)
    for i in range(6):
        for key in names[i % 3:] + names[:i % 3]:
            mat, v = matrices[key], inputs[key]
            mat @ v
            repeats = max(20, min(500, 3_000_000 // mat.nnz))
            start = perf_counter()
            for _ in range(repeats):
                mat @ v
            times[key].append((perf_counter() - start) * 1e6 / repeats)
    for key in times:
        result[key]["matvec_median_us"] = round(median(times[key]), 3)
    return result


def run_width(width, verify):
    receipt_file = OUT/f"width-{width}.json"
    vectors_file = OUT/f"width-{width}-witnesses.npz"
    t0 = perf_counter()
    parts = matrices(width, 6)
    t1 = perf_counter()
    full = matrices(width + 2, 1)
    t2 = perf_counter()
    models = {"base": {1: parts[1]},
              **{f"short{cap}": {d: m for d, m in parts.items() if d <= cap}
                 for cap in range(3, 7)}, "full": full}
    if verify:
        receipt = json.loads(receipt_file.read_text())
        assert sha256(vectors_file.read_bytes()).hexdigest() == receipt["witness_sha256"]
        assert sha256(SOURCE.read_bytes()).hexdigest() == receipt["source_sha256"]
        with np.load(vectors_file, allow_pickle=False) as archive:
            for key, subset in models.items():
                proof = receipt["models"][key]
                margin = check_enclosure(subset, archive[key], proof["site_numerators"])
                assert margin == proof["strict_scaled_slack"]
                a, b = proof["site_numerators"]
                assert proof["order_interval"] == [str(Fraction(t, DEN)**2) for t in (a,b)]
                if key.startswith("short"):
                    gains = gain_fraction(receipt["models"]["base"], proof, receipt["models"]["full"])
                    assert list(map(str, gains)) == proof["gain_fraction_interval"]
        print("verified", width, len(models), flush=True)
        return
    results, vectors = {}, {}
    for key, subset in models.items():
        vectors[key], results[key] = certify(subset)
    for key in results:
        if key.startswith("short"):
            gains = gain_fraction(results["base"], results[key], results["full"])
            results[key]["gain_fraction_interval"] = list(map(str, gains))
            results[key]["gain_percent_approx"] = round(50 * float(sum(gains)), 6)
    np.savez_compressed(vectors_file, **vectors)
    receipt = {
        "width": width, "outside_width": width + 2, "models": results,
        "edge_counts_by_length": {d: {"nonzero": m.nnz, "words": int(m.sum())}
                                  for d, m in parts.items()},
        "exact_coefficient_csr_bytes": {str(d): storage(m) for d, m in parts.items()},
        "single_sweep": sparse_bench(parts[1], parts, full[1], math.sqrt(results["short6"]["order_approx"])),
        "build_seconds": {"inner_and_short_matrices": round(t1-t0, 4), "full_matrix": round(t2-t1, 4)},
        "total_process_seconds": round(perf_counter()-t0, 4),
        "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "witness_sha256": sha256(vectors_file.read_bytes()).hexdigest(),
        "environment": {"python": platform.python_version(), "scipy": scipy.__version__, "numpy": np.__version__,
                        "platform": sys.platform, "numeric_threads": 1},
    }
    receipt_file.write_text(json.dumps(receipt, indent=2) + "\n")
    print(width, {k: x["order_approx"] for k,x in results.items()}, results["short6"]["gain_percent_approx"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=[8,10,12,14])
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for w in args.widths:
        if w < 4 or w % 2:
            parser.error("supply even widths >= 4")
        run_width(w, args.verify)
