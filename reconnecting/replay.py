"""Construct both large certificates, or independently audit existing vectors."""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
EXT = HERE.parent / "extension"
OLD = "23fdd8c621fedbb2ef40556a1bf484758149ed28afebee0aed0a0675123ef538"
SEED = "d7bd9f5549906c7bf5a7e1d591cd083e4c33a511d51369802c66cfb8c03802a2"


def matched(path, expected):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4<<20), b""):
            h.update(block)
    if h.hexdigest()!=expected:
        raise ValueError(f"Witness checksum differs: {path.name}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=8)
    parser.add_argument("--from-v11",type=Path,help="A verified v1.1 round-440 witness")
    parser.add_argument("--audit-vectors",type=Path,nargs=3,metavar=("ROUND500","MIXED","UNIT"))
    parser.add_argument("--windows-only",action="store_true")
    args=parser.parse_args()
    if args.threads<1:
        raise ValueError("Threads must be positive")
    work=args.work.resolve()
    work.mkdir(exist_ok=True,parents=True)
    env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),OMP_DYNAMIC="FALSE")

    def run(*cmd):
        return subprocess.check_output([str(item) for item in cmd],env=env,text=True)

    def compile(name,src,standard="c++17"):
        exe=work/name
        run("g++","-O3","-std="+standard,"-fopenmp","-I",EXT/"research/review73/growth",src,"-o",exe)
        return exe

    if args.windows_only:
        generator=compile("window-build",HERE/"native/window/certify.cpp","c++20")
        audit=compile("window-audit",HERE/"native/window/audit.cpp","c++20")
        for k,p,q,hi in [(12,159,100,38785544598),(13,36,25,37707676858)]:
            prefix=work/f"window-{k}"
            record=json.loads(run(generator,k,p,q,110,10**10,prefix,"--fold"))
            rows=json.loads(run(audit,prefix,k,p,q,hi,10**10,"--fold"))
            assert record["upper"]==hi and rows["all_rows_pass"] and rows["states"]==record["states"]
            print(k,record["upper"],rows["states"])
        return

    if args.audit_vectors:
        seed, mixed, unit = (f.resolve() for f in args.audit_vectors)
    else:
        v11=args.from_v11
        if v11 is None:
            base=work/"v11"
            run(sys.executable,"-B",EXT/"research/review74/boundary/replay.py","--full","--work",base,"--threads",args.threads)
            run(sys.executable,"-B",EXT/"research/review75/growth/build.py",work/"gap")
            v11=work/"round440.bin"
            run(work/"gap/compute",base/"short-exits-32.bin",80,v11)
        matched(v11,OLD)
        seed,mixed,unit=(work/name for name in ("round500.bin","mixed560.bin","unit540.bin"))
        for name,in_file,rounds,out in [("nested",v11,60,seed),("mixed",seed,60,mixed),("unit",seed,40,unit)]:
            if out.exists():
                matched(out,SEED if name=="nested" else json.loads((HERE/"certificates"/(name+".json")).read_bytes())["output_sha256"])
                continue
            generator=compile(name+"-build",HERE/"native"/name/"compute.cpp")
            print(run(generator,in_file,rounds,out),end="")
    matched(seed,SEED)
    env["LRWR_GROWTH77_PREVIOUS"]=str(seed)
    for name,vector in [("mixed",mixed),("unit",unit)]:
        reference=json.loads((HERE/"certificates"/(name+".json")).read_bytes())
        matched(vector,reference["output_sha256"])
        audit=compile(name+"-audit",HERE/"native"/name/"audit.cpp")
        log=run(audit,vector)
        extrema=re.search(r"min (\d+/\d+) max (\d+/\d+)",log)
        assert (Fraction(extrema[1]),Fraction(extrema[2]))==(Fraction(reference["minimum"]),Fraction(reference["maximum"]))
        assert "odd 386058577 even 424317073" in log
        print(name,log)


if __name__=="__main__":
    main()
