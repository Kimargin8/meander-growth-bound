"""Check transcripts or recompute a positive-vector certificate from source."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
import os
from pathlib import Path
import re
import subprocess
import sys
from streamed_path import scan

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "output/research74/boundary"
ORIGINAL = ROOT / "research/review73/growth/implicit_bound.cpp"
CERTIFIED = Fraction(5660040998411745, 483083866680971)
OLD_UPPER = Fraction(15264493432827, 1305926315603)


def last_minimum(text):
    found = re.findall(r"\bmin (\d+/\d+)", text)
    if not found:
        raise ValueError("No minimum recorded")
    return Fraction(found[-1])


def checked_run(arguments, threads):
    return subprocess.run([str(x) for x in arguments], check=True, text=True,
                          stdout=subprocess.PIPE,
                          env={**os.environ,"OMP_NUM_THREADS":str(threads)}).stdout


def file_digest(path):
    check = sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(8*1024**2), b""):
            check.update(part)
    return check.hexdigest()


def evaluate_stored():
    recorded = (DATA/"w32.txt").read_text()
    audited = (DATA/"w32-audit.txt").read_text()
    assert (DATA/"replayed.txt").read_text().splitlines() == recorded.splitlines()
    assert (DATA/"replayed-audit.txt").read_text().splitlines() == audited.splitlines()
    assert "odd 386058577 even 424317073" in recorded and "odd 386058577 even 424317073" in audited
    assert last_minimum(recorded) == last_minimum(audited) == CERTIFIED
    assert Fraction("11.7164769697") < CERTIFIED < 12
    assert OLD_UPPER < CERTIFIED
    assert Fraction(1,12)+Fraction(3,12**2)+Fraction(9,12**3) == Fraction(7,64)
    import json
    from adjacent_weights import check as check_adjacent
    assert check_adjacent() == json.loads((DATA/"adjacent-weights.json").read_bytes())
    previous = json.loads((DATA/"streamed-path.json").read_text())["rows"]
    for entry in previous:
        new = scan(entry["width"],7,2)
        assert all(new[k] == value for k,value in entry.items() if not k.endswith("seconds"))


def native(width, directory, threads):
    directory.mkdir(parents=True, exist_ok=True)
    subprocess.run([sys.executable, str(HERE/"build_weighted.py"), str(directory)], check=True)
    if sys.byteorder != "little":
        raise ValueError("The published starting-vector format is little-endian")
    subprocess.run(["g++","-O3","-std=c++17","-fopenmp",str(ORIGINAL),
                    "-o",str(directory/"original")],check=True)
    base = directory/f"baseline-{width}.bin"
    final = directory/f"short-exits-{width}.bin"
    rounds, extra = (300,60) if width==32 else (100,10)
    # Reuse a baseline only when the known full-size digest has been supplied.
    digest_path = DATA/"checkpoint-sha256.txt"
    digests = dict(line.split()[::-1] for line in digest_path.read_text().splitlines())
    reuse = width == 32 and base.is_file() and file_digest(base) == digests["width32.local.checkpoint"]
    if not reuse:
        checked_run([directory/"original", width, rounds, base], threads)
        if width == 32:
            assert file_digest(base) == digests["width32.local.checkpoint"]
    text = checked_run([directory/"weighted", base, extra, final],threads)
    audit = checked_run([directory/"audit",final],threads)
    value = last_minimum(text)
    assert value == last_minimum(audit)
    if width == 32:
        assert value == CERTIFIED
        assert file_digest(final) == digests["width32-round360.bin"]
    print(f"Width {width}, two independent sweeps agree: {value}")


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--small", action="store_true", help="Build and run both native programs at width 14")
    group.add_argument("--full", action="store_true", help="Regenerate the width-32 certificate (Linux/GCC)")
    parser.add_argument("--work",type=Path, help="Directory for native binaries and two private checkpoints")
    parser.add_argument("--threads",type=int,default=8)
    args=parser.parse_args()
    evaluate_stored()
    print("Stored arithmetic, old-strip comparison, and all-duration outer-layer tests agree.")
    if args.small or args.full:
        if args.work is None or args.threads < 1:
            parser.error("Supply a writable --work directory and positive --threads")
        native(32 if args.full else 14,args.work.resolve(),args.threads)
