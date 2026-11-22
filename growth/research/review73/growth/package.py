"""Build and independently exercise the private, portable review archive."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT/"output/research73"
GROWTH = "research/review73/growth/"
LOGS = "output/research73/growth/"
RETURNS = "research/review73/growth_return/"
RETURN_LOGS = "output/research73/growth_return/"
SHARED = "research/review57/strip/confined.cpp"
KNOWN = "references/oeis/2026-11-15/b005315.txt"
REPLAY = LOGS+"clean-replay.json"
PAPER_SOURCE = "research/review73/star/growth_constant.tex"
PAPER_PDF = "output/research73/paper/growth_constant.pdf"
ARCHIVE = OUT/"paper/growth-constant-review.zip"
CORE_TEST_INPUTS = [GROWTH+f for f in ("verify.py","strip_bound.cpp","folded_bound.cpp","implicit_bound.cpp")]
CORE_TEST_INPUTS += [SHARED,LOGS+"w26.txt",LOGS+"f28.txt"]
if (HERE/"check_packet.py").is_file():
    from check_packet import driver_fingerprint, matching_replay_input
else:
    sys.path.insert(0,str(ROOT))
    from verify_package import driver_fingerprint, matching_replay_input


def digest(blob):
    return sha256(blob).hexdigest()


def members(*,for_replay=False):
    # All scientific programs stay in their original hierarchy, including
    # their relative includes. There are no binaries or generated checkpoints.
    items = {rel:(ROOT/rel).read_bytes() for rel in CORE_TEST_INPUTS}
    if for_replay:
        return items
    files = [KNOWN,LOGS+"growth-bounds.json",REPLAY,GROWTH+"README.md",GROWTH+"package.py",
             GROWTH+"check_strip_edges.cpp",LOGS+"edges32.txt",GROWTH+"finite_bound.cpp",LOGS+"finite32.txt",
             GROWTH+"verify_vector.cpp",LOGS+"vector32-final.txt",
             "research/review73/star/height_envelope.py","research/review73/star/approximation.py",
             "research/review73/star/strip_error.py",
             "output/research73/star/approximation.json",GROWTH+"accuracy.py",LOGS+"accuracy.json",
             RETURNS+"return.cpp",RETURNS+"verify.py",RETURNS+"reading.md",
             RETURN_LOGS+"verified.json",RETURN_LOGS+"return26-time.txt",
             PAPER_SOURCE,PAPER_PDF]
    regular = json.loads((ROOT/LOGS/"growth-bounds.json").read_bytes())
    files.extend(regular["sources_and_transcripts"])
    return_cert = json.loads((ROOT/RETURN_LOGS/"verified.json").read_bytes())
    files.extend(RETURN_LOGS+name for name in return_cert["sha256"] if name.endswith(".txt"))
    items.update({rel:(ROOT/rel).read_bytes() for rel in set(files)})
    readme = HERE/"review-readme.md"
    verifier = HERE/"check_packet.py"
    items["README.md"] = (readme if readme.is_file() else ROOT/"README.md").read_bytes()
    items["verify_package.py"] = (verifier if verifier.is_file() else ROOT/"verify_package.py").read_bytes()
    return items


def write_archive(target,mapping):
    mapping = dict(mapping)
    mapping["SHA256SUMS"] = "".join(
        f"{digest(blob)}  {name}\n" for name,blob in sorted(mapping.items())
    ).encode("ascii")
    with ZipFile(target,"w",compression=ZIP_DEFLATED,compresslevel=9) as archive:
        for name,blob in sorted(mapping.items()):
            info = ZipInfo(name,(1980,1,1,0,0,0))
            info.create_system=3
            info.external_attr=(0o100644 << 16)
            info.compress_type=ZIP_DEFLATED
            archive.writestr(info,blob,compresslevel=9)


def extracted(path,target):
    with ZipFile(path) as archive:
        if archive.testzip():
            raise ValueError("invalid ZIP member checksum")
        for name in archive.namelist():
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("unexpected archive member")
        archive.extractall(target)


def replay(cxx, threads):
    inputs=members(for_replay=True)
    with TemporaryDirectory() as tmp:
        folder=Path(tmp)
        path=folder/"sources.zip"
        root=folder/"extracted"
        write_archive(path,inputs)
        extracted(path,root)
        # Only these files have been placed in the fresh tree. The published
        # verifier compiles all three programs into its own temporary directory.
        command=[sys.executable,"research/review73/growth/verify.py","--small-native",
                 "--replay","w26","--replay","f28","--threads",str(threads),"--cxx",cxx]
        completed=subprocess.run(command,cwd=root,capture_output=True,text=True,check=True)
        if not all(s in completed.stdout for s in ("All three transfers","w26","f28")):
            raise ValueError("incomplete successful replay output")
        receipt=dict(
            native_semantic_reproduction=["w26: full physical graph, 250 iterations",
                                          "f28: reflected graph, 300 iterations"],
            direct_small_reproduction="three freshly compiled transfers, all 7 even widths 2..14, orders 1..8",
            method="extract a new source ZIP into an empty temporary directory and compile there; compare counts, return tests and every rational enclosure",
            platform=platform.system(),python=platform.python_version(),
            compiler=subprocess.check_output([cxx,"--version"],text=True).splitlines()[0],
            threads=threads,
            command=["python3",*command[1:]],
            source_and_expected_log_sha256={k:digest(v) for k,v in sorted(inputs.items())},
            native_replay_driver_sections=driver_fingerprint(inputs[GROWTH+"verify.py"]),
            output=completed.stdout.splitlines())
        dest=ROOT/REPLAY
        tmp_dest=dest.with_suffix(".pending")
        tmp_dest.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf8")
        tmp_dest.replace(dest)
        print(completed.stdout,end="")


def build():
    mapping=members()
    manuscript=mapping[PAPER_SOURCE].decode("utf8")
    if any(token in manuscript for token in ("GROWTHDECIMAL", "GROWTHNUM", "GROWTHDEN")):
        raise ValueError("The manuscript is still awaiting its final exact numbers")
    checked=json.loads(mapping[REPLAY])
    if any(not matching_replay_input(k,mapping[k],checked) for k in checked["source_and_expected_log_sha256"]):
        raise ValueError("The signed native inputs have changed; rerun --replay in a clean directory.")
    ARCHIVE.parent.mkdir(exist_ok=True,parents=True)
    pending=ARCHIVE.with_suffix(".pending")
    write_archive(pending,mapping)
    with TemporaryDirectory() as tmp:
        root=Path(tmp)
        extracted(pending,root)
        subprocess.run([sys.executable,"verify_package.py"],cwd=root,check=True)
    pending.replace(ARCHIVE)
    print(ARCHIVE.name, len(mapping)+1,"files", ARCHIVE.stat().st_size,"bytes",digest(ARCHIVE.read_bytes()))


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Run Python without -O")
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replay",action="store_true",help="compile in isolation and save the replay result")
    ap.add_argument("--threads",type=int,default=6)
    ap.add_argument("--cxx",default="g++")
    args=ap.parse_args()
    if args.threads<1:
        ap.error("threads must be positive")
    if args.replay:replay(args.cxx,args.threads)
    else:build()
