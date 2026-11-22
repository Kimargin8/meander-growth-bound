"""Audit the endpoint of a one-gap detour, and certify three shift-only models."""

from collections import Counter
from fractions import Fraction
import json
from math import comb, sqrt
from hashlib import sha256
from pathlib import Path
import sys

from scan import (DEN, OUT, approximate, certify, check_enclosure, gain_fraction, matchings, matrices)
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import coo_matrix


def literal_splice_after_open(p, gap, new_gap):
    """Glue one cap on labelled ends and follow each surviving arc by its endpoint."""
    top = len(p)
    labels = list(range(gap)) + [top, top+1] + list(range(gap, top))
    partners = dict(enumerate(p))
    partners.update({top:top+1, top+1:top})
    left, right = labels[new_gap-1:new_gap+1]
    if partners[left] == right:
        return None
    surviving = labels[:new_gap-1] + labels[new_gap+1:]
    renumber = {label:i for i,label in enumerate(surviving)}
    endpoints = []
    for end in surviving:
        path = partners[end]
        if path == left:
            path = partners[right]
        elif path == right:
            path = partners[left]
        endpoints.append(renumber[path])
    return tuple(endpoints), new_gap-1


def boundary_shift(width):
    n = sum((2*k+1)*len(matchings(k)) for k in range(1,width//2+1))
    addresses = list(matchings(width//2))
    first = n - (width+1)*len(addresses)
    rows, cols = [], []
    minima = {}
    for gap in range(width+1):
        counts = Counter({gap+1: 1})
        for steps in range(1,6):
            move = Counter()
            for loc, ways in counts.items():
                if loc > 0:
                    move[loc-1] += ways
                if loc < width+2:
                    move[loc+1] += ways
            counts = move
            if steps % 2:
                for delta in (-1,1):
                    if 0 <= gap + delta <= width:
                        count = counts[gap+1+delta]
                        minima[steps] = min(count,minima.get(steps,count))
    assert minima == {1:1,3:3,5:9}
    for i,p in enumerate(addresses):
        offset = first+i*(width+1)
        for gap in range(width+1):
            for delta in (-1,1):
                if not 0 <= gap+delta <= width:
                    continue
                assert literal_splice_after_open(p,gap,gap+1+delta) == (p,gap+delta)
                rows.append(offset+gap)
                cols.append(offset+gap+delta)
    assert len(rows) == 2*width*comb(width,width//2)//(width//2+1)
    S = coo_matrix((np.ones(len(rows),dtype=np.int32),(rows,cols)),shape=(n,n)).tocsr()
    return S, minima


def exact_margins(T,S,v,mode,ends):
    Tv, Sv = T @ v, S @ v
    assert 4*int(max(v)) < 2**63
    values = []
    for i,a in enumerate(ends):
        if mode == "dyadic":
            alpha, beta, gamma = DEN*64, DEN*7, 64*a
        else:
            res = a*a - 3*DEN*DEN
            assert res > 0
            alpha, beta, gamma = DEN*res, DEN**3, a*res
        residuals = [alpha*int(t) + beta*int(s) - gamma*int(c) for t,s,c in zip(Tv,Sv,v)]
        values.append(min(residuals) if i==0 else -max(residuals))
    assert min(values)>0
    return [str(i) for i in values]


def smooth_model(T,S,mode):
    def eig(r):
        coeff = 7/64 if mode == "dyadic" else 1/(r*r-3)
        return approximate({1:T.astype(float) + S*coeff})
    if mode=="dyadic":
        r,v = eig(3)
    else:
        r = brentq(lambda site:eig(site)[0]-site, sqrt(3)+.1, 4.1, xtol=3e-14)
        _,v = eig(r)
    assert mode != "dyadic" or r*r < 12
    center = int(r*DEN)
    ends = [center-1,center+2]
    slack = exact_margins(T,S,v,mode,ends)
    return v, {"site_numerators":ends, "site_denominator":DEN,"strict_scaled_slack":slack,
               "order_approx":r*r, "order_interval":[str(Fraction(a,DEN)**2) for a in ends]}


def main(verify):
    report_file = OUT/"shift-only.json"
    arrays_file = OUT/"shift-only-witnesses.npz"
    if verify:
        report = json.loads(report_file.read_text())
        if sha256(arrays_file.read_bytes()).hexdigest() != report["witness_sha256"]:
            raise ValueError("unexpected witness archive")
        report_data = report["rows"]
    else:
        report_data = {}
    vectors = {}
    for width in (2,4,6,8,10,12,14):
        parts = matrices(width, 3)
        T, E3 = parts[1],parts[3]
        S, minima = boundary_shift(width)
        assert (E3!=S).nnz==0
        if width < 8:
            continue
        original = json.loads((OUT/f"width-{width}.json").read_text())["models"]
        rows = {"labelled_top_shifts": S.nnz, "walk_minima":minima, "models":{}}
        truncated = {1:T, 3:S, 5:S*3, 7:S*9}
        for mode in ("dyadic","truncated","p5"):
            key = f"{width}-{mode}"
            if verify:
                record = report_data[str(width)]["models"][mode]
                with np.load(arrays_file,allow_pickle=False) as archive:
                    witness = archive[key]
                margins = (check_enclosure(truncated,witness,record["site_numerators"]) if mode=="truncated" else
                           exact_margins(T,S,witness,mode,record["site_numerators"]))
                assert margins == record["strict_scaled_slack"]
                assert record["order_interval"] == [str(Fraction(a,DEN)**2) for a in record["site_numerators"]]
            elif mode == "truncated":
                vectors[key], record = certify(truncated)
            else:
                vectors[key], record = smooth_model(T,S,mode)
            gains = gain_fraction(original["base"],record,original["full"])
            if verify:
                assert record["gain_fraction_interval"] == list(map(str,gains))
            record["gain_fraction_interval"] = list(map(str,gains))
            record["gain_percent_approx"] = round(50*float(sum(gains)),6)
            rows["models"][mode] = record
        if not verify:
            report_data[str(width)] = rows
        print(width,{k:(round(v["order_approx"],8),v["gain_percent_approx"])
                     for k,v in rows["models"].items()})
    if not verify:
        np.savez_compressed(arrays_file,**vectors)
        report = {"site_denominator":DEN, "rows":report_data,
                  "witness_sha256":sha256(arrays_file.read_bytes()).hexdigest(),
                  "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}
        report_file.write_text(json.dumps(report,indent=2)+"\n")
    else:
        assert report["source_sha256"] == sha256(Path(__file__).read_bytes()).hexdigest()


if __name__=="__main__":
    main("--verify" in sys.argv)
