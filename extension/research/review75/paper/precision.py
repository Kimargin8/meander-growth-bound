"""Select a certified width and optionally return the count in that strip."""
from argparse import ArgumentParser
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/"output/research75/growth/verified.json"
sys.path.insert(0,str(ROOT/"research/review73/star"))
sys.path.insert(0,str(ROOT/"research/review75/growth"))
from verify import OLD, OLD_C, finite_width


def choose(n,eps):
    if n<1 or not 0<eps<1:
        raise ValueError("Expected a positive order and a relative accuracy in (0,1)")
    cert=json.loads(DATA.read_bytes())
    new=Fraction(cert["x"])
    new_c=Fraction(cert["constant_for_all_orders"])
    if n<=5:
        return 2*n
    return min(finite_width(n,eps,OLD_C,OLD),finite_width(n,eps,new_c,new),
               finite_width(n,eps,Fraction(1,2**21),Fraction(292,25)))


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("order",type=int)
    parser.add_argument("epsilon",type=Fraction)
    parser.add_argument("--count",action="store_true",help="Also execute the ordinary confined count; cost grows exponentially with width")
    opts=parser.parse_args()
    w=choose(opts.order,opts.epsilon)
    answer={"order":opts.order,"relative_error_at_most":str(opts.epsilon),"width":w}
    if opts.count:
        from strip_error import connected_strip
        answer["lower_count"]=connected_strip(w,(opts.order,))[opts.order]
    print(json.dumps(answer))

