"""Find the positive root of the first-return polynomial on the original addresses."""

import math
from scan import DEN, check_enclosure
import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import eigs


def certify_compact(parts):
    # r is measured per site. The operator below equals r B(1/r).
    state_vector = np.ones(parts[1].shape[0])
    calls = 0

    def mismatch(r):
        nonlocal state_vector, calls
        matrix = parts[1].astype(np.float64)
        for length, coefficient in parts.items():
            if length > 1:
                matrix += coefficient * (r**(1 - length))
        eig, vec = eigs(matrix, k=1, which="LR", v0=state_vector, ncv=18,
                        tol=2e-13, maxiter=4000)
        state_vector = vec[:, 0].real.copy()
        if state_vector[0] < 0:
            state_vector *= -1
        assert np.min(state_vector) > 0
        calls += 1
        return eig[0].real - r

    # Four choices at each site, while the strip at W >= 8 has multiple returns.
    r = brentq(mismatch, 1.5, 4.01, xtol=5e-14)
    mismatch(r)
    witness = np.rint(state_vector / state_vector.max() * (1 << 48)).astype(np.int64)
    center = math.floor(r*DEN)
    site_interval = [center-1, center+2]
    slack = check_enclosure(parts, witness, site_interval)
    return witness, {"site_numerators": site_interval, "strict_scaled_slack": slack,
                     "parameter_evaluations": calls}
