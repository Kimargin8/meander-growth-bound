"""Build the two integer transfers using independently computed dyadic edge weights."""
from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIG = ROOT/"research/review73/growth"
TARGET = Fraction(1173, 100)
DEN = 4096


def determinant_edges(m, x):
    d = [1, 1]
    for j in range(2,m+1):
        d.append((x if j % 2 == 0 else 1)*d[-1]-d[-2])
    return [Fraction(d[i]*d[m-i-2]*(x if i%2 and (m-i-2)%2 else 1),d[m]) for i in range(m-1)]


@lru_cache(None)
def walk_weight(width, start, target):
    m = width+3
    v = [0]*m
    v[start] = 1
    summed = Fraction(0)
    for j in range(1,100):
        v = [a+b for a,b in zip([0]+v[:-1],v[1:]+[0])]
        if j % 2:
            h = (j+1)//2
            summed += Fraction(v[target],TARGET**h)
            tail = 2/TARGET*(4/TARGET)**h/(1-4/TARGET)
            lo, hi = int(DEN*summed), int(DEN*(summed+tail))
            if lo == hi:
                return lo
    raise ArithmeticError("Walk enclosure did not isolate a dyadic weight")


def tables():
    cofactor, walk, restricted = [[[]] for _ in range(3)]
    for w in range(2,33,2):
        m = w+3
        exact = determinant_edges(m,TARGET)[1:-1]
        first = [int(DEN*x) for x in exact]
        second = [walk_weight(w,i,i+1) for i in range(1,w+1)]
        # Keep all walks within four steps of their eventual target.
        # This geometric restriction bounds the remaining time to close.
        short = []
        for i in range(1,w+1):
            both = []
            for end in (i,i+1):
                left, right = max(0,end-4), min(m-1,end+4)
                bounded=int(DEN*determinant_edges(right-left+1,TARGET)[i-left])
                assert bounded==walk_weight(right-left-2, i-left, i-left+1)
                both.append(bounded)
            short.append(min(both))
        assert first == first[::-1] == second and max(first)<512
        cofactor.append(first)
        walk.append(second)
        restricted.append(short)
    return cofactor,walk,restricted


def header(rows):
    return ("constexpr unsigned long long gap_weight[17][32]={\n" +
            ",\n".join("{"+",".join(map(str,row))+"}" for row in rows)+"};\n")


def replace_once(code, old, new):
    if code.count(old) != 1:
        raise RuntimeError("Missing or ambiguous insertion point: "+old)
    return code.replace(old,new)


def build(directory):
    direct, series, restricted = tables()
    directory.mkdir(parents=True,exist_ok=True)
    for filename, program, rows, digest in (
        ("implicit_bound.cpp","compute",direct,"f8c885bbe70bbb83ac0aa2a9cd797c3f0c0fe1cc7689fd8eff2eb82b0baac479"),
        ("verify_vector.cpp","audit",series,"b0656177a8b8460cd6a3348a2506cbc90b98c74a8aafe8a1bea6158c8d2e980a")):
        code = (ORIG/filename).read_text()
        assert sha256(code.encode()).hexdigest()==digest
        if program=="compute":
            code = replace_once(code,"                        z+=v[side]>>readshift;",
                """                        U right=v[side]>>readshift;
                        if(k==cap){
                            U left_bonus=u ? z*gap_weight[k][u-1] : 0;
                            U right_bonus=right*gap_weight[k][u];
                            z += (left_bonus+right_bonus)/4096;
                        }
                        z+=right;""")
            code = code[:code.index("int main(int argc,")] + (ROOT/"research/review74/boundary/weighted_main.inc").read_text()
        else:
            start=code.index("                        if (u) z+=src")
            end=code.index("                        if (k<cap)",start)
            code = code[:start]+"""                        U left = u ? src[l.index(1-p,r,u-1)] : 0;
                        U right = u<k ? src[l.index(1-p,r,u+1)]
                                      : src[l.index(1-p,l.rank(reflected),u-1)];
                        z = left+right;
                        if (k==cap){
                            Wide bonus = Wide(right)*gap_weight[k][u];
                            if(u) bonus += Wide(left)*gap_weight[k][u-1];
                            z += U(bonus>>12);
                        }
"""+code[end:]
            code = replace_once(code,'if (argc<2 || argc>3 || (argc==3 && std::string(argv[2])!="--shift"))',
                                    'if (argc!=2)')
        generated = directory/(program+".cpp")
        generated.write_text(header(rows)+code)
        subprocess.run(["g++","-O3","-fopenmp","-std=c++17","-I",str(ORIG),str(generated),
                        "-o",str(directory/program)],check=True)
    print(json.dumps({"squared_parameter":str(TARGET),"denominator":DEN,
                      "width32_all":direct[16], "width32_radius4":restricted[16],
                      "all_tested_widths":list(range(2,33,2))}))


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    build(p.parse_args().directory)
