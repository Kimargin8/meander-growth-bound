"""Literal checks for sublinear excess, subdivision and codes of loop-free halves.

Words run from the cut towards the fresh vertices. The initial height is
determined by the number of down-steps minus up-steps on that bank.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from math import comb
from pathlib import Path
import json


@cache
def words(h,t):
    if h<0 or h>t or (h-t)%2:return ()
    if t==0:return ("",)
    return tuple(c+w for c,change in (("U",1),("D",-1)) for w in words(h+change,t-1))


@cache
def arcs(word):
    h=len(word)-2*word.count("U")
    if h<0:return None
    stack=list(range(-h,0));mates=[None]*len(word)
    for i,c in enumerate(word):
        if c=="U":stack.append(i)
        elif stack:
            j=stack.pop();mates[i]=j
            if j>=0:mates[j]=i
        else:return None
    return tuple(mates) if not stack else None


def components(up,low):
    banks=arcs(up),arcs(low)
    if None in banks:return None
    seen=set();out=[]
    for start in range(len(up)):
        if start in seen:continue
        todo=[start];seen.add(start);sites=[];stubs=[]
        while todo:
            i=todo.pop();sites.append(i)
            for side,bank in enumerate(banks):
                j=bank[i]
                if j<0:stubs.append((side,-j))
                elif j not in seen:seen.add(j);todo.append(j)
        if not stubs:return None
        assert len(stubs)==2
        least=1 if stubs[0][0]!=stubs[1][0] else 2
        extra=len(sites)-least
        assert extra>=0 and extra%2==0
        out.append((tuple(sorted(stubs)),tuple(sorted(sites)),extra//2))
    return out


def signature(parts):
    return tuple(sorted(stubs for stubs,_,_ in parts))


def full_connected(up,low):
    u,v=arcs(up),arcs(low)
    assert u and v and all(i>=0 for i in (*u,*v))
    seen={0};stack=[0]
    while stack:
        i=stack.pop()
        for j in (u[i],v[i]):
            if j not in seen:seen.add(j);stack.append(j)
    return len(seen)==len(up)


def subdivide(up,low,i):
    # Upper old edge at the first new vertex; lower old edge at the third.
    return up[:i+1]+"UD"+up[i+1:],low[:i]+"UD"+low[i:]


def prims(length):
    return ["U"+v+"D" for v in words(0,length-2)]


def make_atoms(limit,kind):
    ordinary=[("D","D")]
    for t in range(2,limit+1,2):
        for v in prims(t):
            ordinary.extend([("D"*t,v),(v,"D"*t)])
    specials=[("UDDD","UUDD"),("UUDD","UDDD")]
    if kind in ("eight","forty_four"):
        specials.extend([("UDD","DUD"),("DUD","UDD"),("DUDD","UDUD"),("UDUD","DUDD"),
                         ("DDUD","UUDD"),("UUDD","DDUD")])
    if kind=="forty_four":
        for a,b,t,size in ((1,1,5,8),(2,0,6,14),(0,2,6,14)):
            paths=[(u,v) for u in words(a,t) for v in words(b,t)
                   if len(components(u,v) or [])==1]
            assert len(paths)==size
            specials.extend(paths)
    for u,v in ordinary:
        parts=components(u,v)
        assert parts and sum(d for _,_,d in parts)==0
    for u,v in specials:
        parts=components(u,v)
        assert parts and len(parts)==1 and parts[0][2]>=1
    if kind=="ordinary":specials=[]
    return [*ordinary,*((u,v) for u,v in specials if len(u)<=limit)]


def verify_code(limit,kind):
    atoms=make_atoms(limit,kind)
    trie={}
    for top,low in atoms:
        node=trie
        for pair in zip(top,low):
            assert "end" not in node
            node=node.setdefault(pair,{})
        assert not node
        node["end"]=len(top)
    generated=[{("","")}]+[set() for _ in range(limit)]
    scalars=[1]
    weights=Counter(len(x) for x,y in atoms)
    for t in range(1,limit+1):
        for u,v in atoms:
            k=len(u)
            if k>t:continue
            for x,y in generated[t-k]:
                w=u+x,v+y
                assert w not in generated[t]
                generated[t].add(w)
        for u,v in generated[t]:
            node=trie;lengths=[]
            for pair in zip(u,v):
                node=node[pair]
                if "end" in node:
                    lengths.append(node["end"]);node=trie
            assert node is trie and sum(lengths)==t
            assert components(u,v) is not None
        scalars.append(sum(f*scalars[t-i] for i,f in weights.items() if i<=t))
        assert scalars[t]==len(generated[t])
    return scalars


def literal_halves(limit,H):
    hist=[]
    through=[]
    squares=[1,2,8,42,262]
    subdivision=Counter();injected=0;deleted=0
    min_words=[[] for _ in range(limit+1)]
    for t in range(limit+1):
        by=Counter()
        one_path=0
        physical=[(h,w) for h in range(t%2,t+1,2) for w in words(h,t)]
        for a,x in physical:
            for b,y in physical:
                parts=components(x,y)
                if a==b==1:
                    assert (parts is not None)==full_connected("U"+x,"U"+y)
                if parts is None:continue
                d=sum(z for _,_,z in parts)
                by[d]+=1
                if d==0:min_words[t].append((x,y))
                if a==b==1:
                    # An old through path becomes a closed meander when a
                    # first fresh vertex with one opening on each side is added.
                    assert len(parts)==1
                    one_path+=1
                removed={i for _,sites,extra in parts if extra for i in sites}
                assert len(removed)<=4*d
                ux="".join(x[i] for i in range(t) if i not in removed)
                vy="".join(y[i] for i in range(t) if i not in removed)
                remainder=components(ux,vy)
                assert remainder is not None and all(z==0 for _,_,z in remainder)
                xx=iter(ux);yy=iter(vy)
                assert "".join(x[i] if i in removed else next(xx) for i in range(t))==x
                assert "".join(y[i] if i in removed else next(yy) for i in range(t))==y
                deleted+=1
                if t<=5:
                    for i in range(t):
                        u,v=subdivide(x,y,i)
                        grafted=components(u,v)
                        assert grafted is not None and signature(grafted)==signature(parts)
                        assert sum(z for _,_,z in grafted)==d+1
                        subdivision[x[i]+y[i]]+=1
        assert by[0]==H[t]
        for d,mass in by.items():
            assert mass<=sum(comb(t,e)*4**e*H[t-e] for e in range(min(t,4*d)+1))
        if t%2:
            assert one_path==squares[t//2]
            through.append([t,(t+1)//2,one_path])
        hist.append([by[i] for i in range(1+max(by))])
    assert set(subdivision)=={"UU","UD","DU","DD"}
    for t in range(1,limit+1):
        for d in range((t-1)//2+1):
            grown=set()
            for x,y in min_words[t-2*d]:
                old=components(x,y)
                for _ in range(d):x,y=subdivide(x,y,0)
                curr=components(x,y)
                assert curr is not None and signature(old)==signature(curr) and sum(z for _,_,z in curr)==d
                assert (x,y) not in grown
                grown.add((x,y));injected+=1
            assert len(grown)==H[t-2*d] and len(grown)<=hist[t][d]
    return dict(excess_histograms=hist,first_point_injections=injected,deleted_and_recovered=deleted,
                subdivision_by_old_letters=dict(sorted(subdivision.items())),through_blocks=through)


def bracket(predicate,places=10):
    # The monotone predicate becomes true immediately to the right of the root.
    scale=10**places;lo=0;hi=4*scale
    assert not predicate(Fraction(lo,scale)) and predicate(Fraction(hi,scale))
    while hi-lo>1:
        mid=(lo+hi)//2
        if predicate(Fraction(mid,scale)):hi=mid
        else:lo=mid
    return dict(lower_numerator=lo,upper_numerator=hi,denominator=scale)


def constants():
    def exceeds(beta,third,fourth,fifth=0,sixth=0):
        if beta<=2:return False
        z=1/beta
        return 1-4*z*z>(z+third*z**3+fourth*z**4+fifth*z**5+sixth*z**6)**2
    A=4*5**12*25**25
    C=24**24
    # Closed meander numbers M10 and M32, OEIS A005315.
    m10=8152860
    m32=154967080111946539652848747626
    assert m10**2>5**19 and A**63<C**63*m32**25
    ratio_power=Fraction(A,C)**63/m32**25
    assert ratio_power<Fraction(931512,10**6)**1575
    separator=Fraction(68,25)
    assert Fraction(A,C)<separator**25 and not exceeds(separator,2,6,8,28)
    return dict(two_extra_atoms=bracket(lambda x:exceeds(x,0,2)),
                eight_extra_atoms=bracket(lambda x:exceeds(x,2,6)),
                forty_four_extra_atoms=bracket(lambda x:exceeds(x,2,6,8,28)),
                entropy_upper_at_one_percent=bracket(lambda x:x**25>Fraction(A,C)),
                fixed_block_M32_lower=bracket(lambda x:x**63>m32),
                fraction_exponential_upper=bracket(lambda x:x**1575>ratio_power),
                M10=m10,M32=m32,exact_inequalities=["M10^2 > 5^19","A^63 < C^63 M32^25",
                                                      "(A/C)^63 / M32^25 < (931512/1000000)^1575"])


def run():
    limit=9
    series={k:verify_code(limit,k) for k in ("ordinary","two","eight","forty_four")}
    return dict(limit=limit,atom_code_counts=series,half_checks=literal_halves(limit,series["ordinary"]),
                constants=constants(),source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=="__main__":
    p=ArgumentParser(description=__doc__);p.add_argument("--write",type=Path);args=p.parse_args()
    data=run()
    if args.write:
        temp=args.write.with_suffix(".pending");temp.write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")
        temp.replace(args.write)
    else:
        script=Path(__file__).resolve()
        path=script.parent/"certificates/excess-growth.json"
        if not path.exists():path=script.parents[3]/"output/research73/star/excess-growth.json"
        assert json.loads(path.read_bytes())==data
    print(data["half_checks"],data["constants"])
