"""Verify an extracted growth-constant review package. Standard Python only."""
import argparse
import ast
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys


def driver_fingerprint(source):
    """Functions and run parameters used by the archived clean native replay."""
    code=source.decode("utf-8").replace("\r\n","\n")
    tree=ast.parse(code)
    names={"compile_one","native_tests","arches","literal","parse"}
    result={node.name:sha256(ast.get_source_segment(code,node).encode("utf8")).hexdigest()
            for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names}
    if set(result)!=names:
        raise ValueError("clean replay functions cannot be identified")
    decl=next(node for node in tree.body if isinstance(node,ast.Assign) and
              any(isinstance(t,ast.Name) and t.id=="RUNS" for t in node.targets))
    runs=ast.literal_eval(decl.value)
    return {"functions":result,"runs":{k:list(runs[k]) for k in ("w26","f28")}}


def matching_replay_input(name,actual,receipt):
    if sha256(actual).hexdigest()==receipt["source_and_expected_log_sha256"][name]:
        return True
    # Later certificate types can add unrelated checks to the same module.
    return (name=="research/review73/growth/verify.py" and
            driver_fingerprint(actual)==receipt.get("native_replay_driver_sections"))


def inventory(root):
    recorded = {}
    for line in (root / "SHA256SUMS").read_text(encoding="ascii").splitlines():
        digest, rel = line.split("  ",1)
        if rel in recorded or Path(rel).is_absolute() or ".." in Path(rel).parts:
            raise ValueError(f"unexpected manifest member {rel}")
        path = root / rel
        actual = sha256(path.read_bytes()).hexdigest()
        if digest != actual:
            raise ValueError(f"SHA-256 mismatch for {rel}")
        recorded[rel] = actual
    return recorded


def exact_fraction(text):
    p,q = map(int,text.split("/"))
    if p <= 0 or q <= 0:
        raise ValueError("All reported ratio numerators and denominators must be positive")
    return Fraction(p,q)


def exact_intersections(root):
    # Combine all reported enclosures for one width, including the
    # quotient and unquotiented implementations at width 26.
    base = root / "output/research73/growth"
    known = json.loads((base/"growth-bounds.json").read_bytes())
    intervals = {}
    for key,trial in known["trials"].items():
        width = trial["width"]
        for line in (base/(key+".txt")).read_text().splitlines():
            parts = line.split()
            if not parts or parts[0] != "round":
                continue
            if (len(parts) != 10 or parts[2] != "min" or parts[6] != "max"
                    or parts[4] not in ("row","node") or parts[8] != parts[4]):
                raise ValueError(f"malformed extremal-ratio record {key}")
            lo,hi = map(exact_fraction,(parts[3],parts[7]))
            if lo > hi:
                raise ValueError(f"reversed ratio in {key}")
            old = intervals.get(width,(lo,hi))
            intervals[width] = max(lo,old[0]),min(hi,old[1])
    if any(lo>hi for lo,hi in intervals.values()):
        raise ValueError("Conflicting bounds for the same width")
    return intervals


def independent_witness(root):
    folder=root/"output/research73/growth"
    original=(folder/"i32.txt").read_text().splitlines()
    other=(folder/"vector32-final.txt").read_text().splitlines()
    stamp=other[0].split()
    if stamp[0:3]!=["width","32","round"] or len(other)!=4:
        raise ValueError("unexpected independent vector transcript")
    m=int(stamp[3])
    observed=other[1].split()
    if observed[::2]!=["min","max"]:
        raise ValueError("missing independently computed bounds")
    rounds=[s.split() for s in original if s.startswith("round ")]
    reference=next(r for r in rounds if int(r[1])==m)
    if m!=max(int(r[1]) for r in rounds)//100*100:
        raise ValueError("independent vector is from a different checkpoint round")
    if [exact_fraction(observed[k]) for k in (1,3)]!=[exact_fraction(reference[k]) for k in (3,7)]:
        raise ValueError("independent width-32 vector does not reproduce both extrema")
    shifts={}
    for row in other[2:]:
        parts=row.split()
        if parts[:1]!=["shift"] or parts[1] in shifts:
            raise ValueError("invalid shifted-vector comparison")
        if parts[2:]==["unavailable"]:
            shifts[parts[1]]=None
        elif len(parts)==6 and parts[2::2]==["min","max"]:
            lo,hi=map(exact_fraction,(parts[3],parts[5]))
            if lo>hi:
                raise ValueError("invalid positive shifted-vector interval")
            shifts[parts[1]]=lo,hi
        else:
            raise ValueError("missing positivity status for a shifted vector")
    if set(shifts)!={"23/2","58/5"}:
        raise ValueError("unexpected set of shifted comparisons")
    certified={}
    for line in (folder/"finite32.txt").read_text().splitlines():
        words=line.split()
        if words[:1]==["candidate"]:
            shift=Fraction(int(words[2]),int(words[1]))
            certified[shift]=exact_fraction(words[4]),exact_fraction(words[6])
    if certified.get(Fraction(0))!=(exact_fraction(reference[3]),exact_fraction(reference[7])):
        raise ValueError("the finite prefactor uses a different ordinary certificate")
    for label, pair in shifts.items():
        if certified.get(Fraction(label))!=pair:
            raise ValueError("independent and finite-prefactor shifted-vector ratios differ")
    return m,shifts


def verify(root):
    if not __debug__:
        raise RuntimeError("Run Python without -O so certificate assertions are active")
    records = inventory(root)
    subprocess.run([sys.executable,"research/review73/growth/verify.py"],cwd=root,check=True)
    subprocess.run([sys.executable,"research/review73/growth_return/verify.py"],cwd=root,check=True)
    subprocess.run([sys.executable,"research/review73/growth/accuracy.py"],cwd=root,check=True)
    small=subprocess.check_output([sys.executable,"research/review73/growth/accuracy.py",
                                   "--order","6","--epsilon","1/100","--count"],cwd=root,text=True)
    if json.loads(small)["lower_count"]!=1828:
        raise ValueError("the executable width-restricted counting example failed")
    intervals = exact_intersections(root)
    checkpoint,shifts=independent_witness(root)
    for comparison in shifts.values():
        if comparison:
            l,u=intervals[32]
            intervals[32]=max(l,comparison[0]),min(u,comparison[1])
    if intervals[32][0]>intervals[32][1]:
        raise ValueError("independent positive witnesses have disjoint enclosures")
    proof = json.loads((root/"output/research73/growth/clean-replay.json").read_bytes())
    if not proof["native_semantic_reproduction"] or not proof["direct_small_reproduction"]:
        raise ValueError("missing clean native or direct small replay")
    for rel, digest in proof["source_and_expected_log_sha256"].items():
        if rel not in records or not matching_replay_input(rel,(root/rel).read_bytes(),proof):
            raise ValueError("replay input has changed: "+rel)
    largest = max(intervals)
    lo,hi = intervals[largest]
    manuscript=(root/"research/review73/star/growth_constant.tex").read_text(encoding="utf8")
    def macro(name):
        match=re.search(r"\\newcommand\{\\"+name+r"\}\{([0-9.]+)\}",manuscript)
        if not match:
            raise ValueError("missing certified manuscript constant "+name)
        return match.group(1)
    if Fraction(int(macro("wnum")),int(macro("wden")))!=lo or Fraction(macro("result"))>=lo:
        raise ValueError("the manuscript growth fraction or strict decimal differs from the exact evidence")
    scale = 10**6
    numerator = lo.numerator*scale//lo.denominator
    print(f"Verified {len(records)} hashes, {len(intervals)} exact width enclosures, and both independent extrema at round {checkpoint}.")
    print(f"Width {largest}: {lo} <= strip growth <= {hi}.")
    print(f"The lower endpoint rounded down to six places is {numerator//scale}.{numerator%scale:06d}.")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root",type=Path,default=Path(__file__).resolve().parent)
    opts=ap.parse_args()
    verify(opts.root.resolve())
