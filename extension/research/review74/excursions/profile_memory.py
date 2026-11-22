"""Measure the full certificate search in fresh serial processes for each variant."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter


def peak():
    if sys.platform == "win32":
        import psutil
        return psutil.Process().memory_info().peak_wset
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)


def worker(width, kind):
    from scan import matrices, certify
    initial = peak()
    start = perf_counter()
    compact = kind == "short6-compact"
    parts = matrices(width+2 if kind=="full" else width, 6 if kind.startswith("short6") else 1)
    built = perf_counter()
    after_build = peak()
    if compact:
        from fixed_parameter import certify_compact
        _, result = certify_compact(parts)
    else:
        _, result = certify(parts)
    end = perf_counter()
    return {"width":width, "model":kind, "start_process_peak_bytes":initial,
            "after_build_peak_bytes":after_build,"total_process_peak_bytes":peak(),
            "build_seconds":round(built-start,6), "certify_seconds":round(end-built,6),
            "site_numerators":result["site_numerators"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width",type=int)
    parser.add_argument("--model",choices=["base","short6","full","short6-compact"])
    args=parser.parse_args()
    if args.model is not None:
        print(json.dumps(worker(args.width,args.model)))
    else:
        root=Path(__file__).resolve().parents[3]
        out=root/"output/research74/excursions"
        results=[]
        names=["base","short6","full","short6-compact"]
        for index,width in enumerate([8,10,12,14]):
            original=json.loads((out/f"width-{width}.json").read_text())["models"]
            for kind in names[index%len(names):]+names[:index%len(names)]:
                value=subprocess.check_output([sys.executable,"-B",__file__,"--model",kind,
                                               "--width",str(width)],text=True)
                measured=json.loads(value)
                original_kind = "short6" if kind=="short6-compact" else kind
                assert measured["site_numerators"]==original[original_kind]["site_numerators"]
                results.append(measured)
        (out/"peak-memory.json").write_text(json.dumps(results,indent=2)+"\n")
        for row in results:
            print(row["width"],row["model"],row["total_process_peak_bytes"],row["after_build_peak_bytes"],
                  round(row["build_seconds"]+row["certify_seconds"],4))
