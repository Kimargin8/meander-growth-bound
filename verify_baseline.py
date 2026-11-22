"""Verify the scientific supplement from an unpacked release, with ordinary Python."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent
COMPANION = ROOT / "companion"
COMPANION_CHECKS = (
    "check_excess_caps.py",
    "check_excess_growth.py",
    "supplement/research/review73/star/check_central_excess.py",
    "height_envelope.py",
    "cut_excess_bound.py",
    "all_boundaries.py",
    "approximation.py",
    "strip_error.py",
    "supplement/research/review73/wholecut/strict_order_kernel.py",
)


def run(cwd, name, *args):
    completed = subprocess.run(
        [sys.executable, "-BS", name, *args], cwd=cwd, capture_output=True, text=True
    )
    if completed.returncode:
        sys.stderr.write(completed.stdout + completed.stderr)
        raise SystemExit(completed.returncode)
    print(f"Passed {name}", flush=True)


def main():
    if not __debug__:
        raise SystemExit("Run Python without -O to enable exact assertions.")
    table = (ROOT / "SHA256SUMS").read_text(encoding="ascii").splitlines()
    paths = set()
    for line in table:
        h, f = line.split("  ", 1)
        path = ROOT / f
        if f in paths or not path.is_relative_to(ROOT) or ".." in Path(f).parts:
            raise ValueError(f"Invalid manifest path {f}")
        if sha256(path.read_bytes()).hexdigest() != h:
            raise ValueError(f"SHA-256 mismatch for {f}")
        paths.add(f)
    print(f"Verified {len(paths)} files.", flush=True)
    run(ROOT / "growth", "verify_package.py")
    with TemporaryDirectory() as directory:
        test = Path(directory) / "minimum.json"
        run(COMPANION, "minimum_caps.py", "--write", str(test))
        expected = COMPANION / "certificates/minimum-caps.json"
        if json.loads(test.read_bytes()) != json.loads(expected.read_bytes()):
            raise ValueError("The shortest-half table differs from the certificate.")
    for name in COMPANION_CHECKS:
        args = ("--audit-root", str(ROOT/"growth")) if name == "strip_error.py" else ()
        run(COMPANION, name, *args)
    print("The growth certificate and all ten companion checks passed.")


if __name__ == "__main__":
    main()
