"""Build two separate physical strip implementations with certified top-slit weights."""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1] / "review73" / "growth"


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    args = parser.parse_args()
    args.build.mkdir(parents=True, exist_ok=True)
    for source, name, anchor, line, expected in (
        ("implicit_bound.cpp", "weighted", "                        if(k<cap){",
         "                        if(k==cap) z += (7*z)/64;",
         "f8c885bbe70bbb83ac0aa2a9cd797c3f0c0fe1cc7689fd8eff2eb82b0baac479"),
        ("verify_vector.cpp", "audit", "                        if (k<cap) {",
         "                        if (k == cap) z += (7*z) >> 6;",
         "b0656177a8b8460cd6a3348a2506cbc90b98c74a8aafe8a1bea6158c8d2e980a"),
    ):
        raw = (BASE / source).read_bytes()
        if sha256(raw.replace(b"\r\n",b"\n")).hexdigest() != expected:
            raise ValueError(f"Unexpected baseline source {source}")
        code = raw.decode().replace("\r\n", "\n")
        if code.count(anchor) != 1:
            raise ValueError(f"Ambiguous top-weight insertion in {source}")
        # Both original loops have only accumulated the two labelled shifts here.
        code = code.replace(anchor, line + "\n" + anchor)
        if name == "weighted":
            code = code[:code.index("int main(int argc,")] + (HERE / "weighted_main.inc").read_text()
        else:
            # The prior shifted-vector audit assumes an exactly linear sweep.
            # Only the ordinary lower-vector comparison applies after rounding.
            code = code.replace("if (argc<2 || argc>3", "if (argc!=2 || argc>3")
        generated = args.build / (name + ".cpp")
        generated.write_text(code)
        subprocess.run(["g++", "-O3", "-fopenmp", "-std=c++17", "-I", str(BASE),
                        str(generated), "-o", str(args.build / name)], check=True)


if __name__ == "__main__":
    main()
