"""Verify this version from an extracted release (standard library, Python 3.10+)."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
EXT=ROOT/"extension"


def run(script,*args):
    return subprocess.check_output([sys.executable,"-B",str(EXT/script),*args],cwd=EXT,text=True)


def main():
    if not __debug__:
        raise RuntimeError("Use ordinary Python with assertions enabled")
    manifests=(ROOT/"SHA256SUMS").read_text(encoding="ascii").splitlines()
    for line in manifests:
        expected,name=line.split("  ",1)
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError(name)
    run("research/review74/boundary/replay.py")
    report=json.loads(run("research/review75/growth/verify.py"))
    assert report["x"]=="12401248195833/1058038773811"
    prefix=json.loads(run("research/review75/proof/check_prefix.py","--x",report["x"],
                          "--root",str(report["root"]),"--maximum",str(report["maximum"])))
    assert prefix==json.loads((EXT/"research/review75/proof/prefix-certificate.json").read_bytes())
    assert json.loads(run("research/review75/paper/precision.py","1000","1/100"))["width"]==624
    assert json.loads(run("research/review75/paper/precision.py","5","1/100","--count"))["lower_count"]==262
    small={(1,"1/2"):1,(2,"1/3"):2,(3,"1/4"):8,(4,"1/5"):42,(6,"1/1000"):1828}
    for (n,eps),count in small.items():
        assert json.loads(run("research/review75/paper/precision.py",str(n),eps,"--count"))["lower_count"]==count
    subprocess.run([sys.executable,"-B",str(ROOT/"reconnecting/certificate.py")],cwd=ROOT,check=True)
    print(f"Verified {len(manifests)} files, boundary weights, growth certificate, finite estimates and small counts.")


if __name__=="__main__":
    main()

