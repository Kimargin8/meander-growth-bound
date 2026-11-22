"""Test minimum-length physical caps and two-sided reconstruction on every cut height.

Only the Python standard library is used. No previous data or algorithm is loaded.
For the reachability condition see Di Francesco, Guitter and Jacobsen,
arXiv:cond-mat/0003008, equation (3.1) and the paragraph following it.
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
def tails(height, length):
    if height < 0 or height > length or (length-height)%2:
        return ()
    if not length:
        return ("",)
    return tuple(c+w for c,d in (("U",1),("D",-1)) for w in tails(height+d,length-1))


def reverse(word):
    return word[::-1].translate(str.maketrans("UD","DU"))


@cache
def pairing(word):
    partner=[-1]*len(word)
    stack=[]
    for i,c in enumerate(word):
        if c=="U":
            stack.append(i)
        else:
            j=stack.pop()
            partner[i]=j
            partner[j]=i
    assert not stack
    return tuple(partner)


def one_cycle(p,q):
    i=0
    for count in range(2,len(p)+1,2):
        i=q[p[i]]
        if not i:
            return count==len(p)
    raise AssertionError("two perfect matchings must close a cycle")


def bank_segment(indices,p):
    lookup={i:j for j,i in enumerate(indices)}
    return "".join("U" if j<lookup[p[i]] else "D" for j,i in enumerate(indices))


def make_cap(a,b,p):
    # Cyclic order at the cut: upper innermost to outermost, then lower
    # outermost to innermost. At new points both stacks pop innermost first.
    up=list(range(a)); low=list(reversed(range(a,a+b)))
    top=[];bottom=[]
    while up or low:
        hi=next((i for i,u in enumerate(up) if p[u]>=a),len(up))
        lo=next((i for i,v in enumerate(low) if p[v]<a),len(low))
        top.extend(("D"*hi,bank_segment(low[:lo],p)))
        bottom.extend((bank_segment(up[:hi],p),"D"*lo))
        up=up[hi:];low=low[lo:]
        if up or low:
            assert up and low and p[up[0]]==low[0]
            top.append("D");bottom.append("D")
            up=up[1:];low=low[1:]
    return "".join(top),"".join(bottom)


def trace_cap(a,b,top,bottom):
    # Independently draw every arc on a new vertex, then walk from each
    # old stub to its partner. A cycle avoiding all stubs spoils a closure.
    k=a+b
    assert len(top)==len(bottom)
    graph=[[] for _ in range(k+len(top))]
    for word,starts in ((top,list(reversed(range(a)))),(bottom,list(range(a,k)))):
        for i,c in enumerate(word,k):
            if c=="U":
                starts.append(i)
            else:
                j=starts.pop()
                graph[i].append(j)
                graph[j].append(i)
        assert not starts
    result=[-1]*k
    internal=0
    for i in range(k):
        if result[i]>=0:continue
        prev=i;cur=graph[i][0]
        while cur>=k:
            internal+=1
            first,second=graph[cur]
            prev,cur=cur,second if first==prev else first
        result[i]=cur
        result[cur]=i
    return tuple(result) if internal==len(top) else None


def inverse(gram):
    n=len(gram)
    rows=[[Fraction(z) for z in r]+[Fraction(i==j) for j in range(n)] for i,r in enumerate(gram)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot]=rows[pivot],rows[j]
        z=rows[j][j]
        rows[j]=[entry/z for entry in rows[j]]
        for i in range(n):
            if i!=j and (z:=rows[i][j]):
                rows[i]=[x-z*y for x,y in zip(rows[i],rows[j])]
    return [r[n:] for r in rows]


def full_rank_minor(matrix):
    work=[[Fraction(x) for x in row] for row in matrix]
    row_ids=list(range(len(work)))
    columns=[]
    for j in range(len(work[0])):
        h=len(columns)
        i=next((i for i in range(h,len(work)) if work[i][j]),None)
        if i is None:continue
        work[h],work[i]=work[i],work[h]
        row_ids[i],row_ids[h]=row_ids[h],row_ids[i]
        for i in range(h+1,len(work)):
            if work[i][j]:
                z=work[i][j]/work[h][j]
                work[i]=[x-z*y for x,y in zip(work[i],work[h])]
        columns.append(j)
    return row_ids[:len(columns)],columns


@cache
def matchings(ends):
    caps=tuple(pairing(w) for w in tails(0,ends))
    assert len(caps)==len(set(caps))==comb(ends,ends//2)//(ends//2+1)
    return caps


@cache
def inverse_gram(ends):
    caps=matchings(ends)
    return inverse([[one_cycle(p,q) for q in caps] for p in caps])


def test_lengths():
    totals=Counter()
    for a in range(7):
        for b in range(7):
            if not a+b or (a-b)%2:continue
            caps=matchings(a+b)
            expected={}
            for p in caps:
                count=sum(j>=a for j in p[:a])
                L=a+b-count
                top,low=make_cap(a,b,p)
                assert len(top)==len(low)==L and trace_cap(a,b,top,low)==p
                expected[p]=L
            minimum={}
            for length in range(max(a,b),max(expected.values())+3,2):
                present=set()
                for u in tails(a,length):
                    for v in tails(b,length):
                        totals["word_pairs_examined"]+=1
                        signature=trace_cap(a,b,u,v)
                        if signature is not None:
                            minimum.setdefault(signature,length)
                            present.add(signature)
                assert present=={p for p,L in expected.items() if L<=length}
            assert minimum==expected
            totals["height_pairs"]+=1
            totals["boundary_matchings"]+=len(expected)
            totals["maximum_cap_length"]=max(totals["maximum_cap_length"],max(expected.values()))
    return dict(totals)


def centroid(word):
    height=last=0;sizes=[]
    for i,c in enumerate(word,1):
        height+=1 if c=="U" else -1
        if not height:
            sizes.append((i-last)//2)
            last=i
    n=len(word)//2;T=sizes[0]
    if 2*T>n+1 or max(sizes)>T:return 0
    first=next((i for i in range(1,len(sizes)) if sizes[i]==T),len(sizes))
    return (1 if 2*T==n+1 else 2)*sum(sizes[:first])


def literal_totals(n,t):
    words=tails(0,2*n)
    items=[(pairing(w),2*w[:t].count("U")-t,centroid(w)) for w in words]
    plain=Counter();weighted=Counter()
    for u,a,_ in items:
        for v,b,weight in items:
            # Flood-fill the ordinary closed graph without choosing an edge
            # colour to follow. The capped counts use a two-colour cycle walk.
            reached=1;todo=[0]
            while todo:
                i=todo.pop()
                for j in (u[i],v[i]):
                    bit=1<<j
                    if not reached&bit:
                        reached|=bit
                        todo.append(j)
            if reached==(1<<(2*n))-1:
                plain[a,b]+=1
                weighted[a,b]+=weight
    return plain,weighted,len(words)**2


def closed_rows(upper,lower,caps,suffix):
    # Count closed meanders with the actual lower half-word kept as a tag.
    # The artificial cap may have a different length for each measurement.
    observed={word:[] for word in lower}
    for x,y in caps:
        if suffix:
            p,q=reverse(x),reverse(y)
            tops=[pairing(p+u) for u in upper]
            bots=[pairing(q+v) for v in lower]
        else:
            tops=[pairing(u+x) for u in upper]
            bots=[pairing(v+y) for v in lower]
        for word,v in zip(lower,bots):
            observed[word].append(sum(one_cycle(u,v) for u in tops))
    return observed


def combine(left,right,inv,whole):
    result={}
    for p,row in left.items():
        true=[sum(u*v for u,v in zip(row,col)) for col in zip(*inv)]
        if whole:assert all(x.denominator==1 and x>=0 for x in true)
        keep=[(i,int(x) if x.denominator==1 else x) for i,x in enumerate(true) if x]
        for u,v in right.items():
            count=sum(x*v[i] for i,x in keep)
            assert count.denominator==1 and count>=0
            result[p,u]=int(count)
    return result


def by_height(a,b,t,R):
    if a+b==0:return [a,b,0,0,0,0,0,0,0,False]
    caps=matchings(a+b)
    lengths=[a+b-sum(j>=a for j in p[:a]) for p in caps]
    left_ids=[i for i,L in enumerate(lengths) if L<=t]
    right_ids=[i for i,L in enumerate(lengths) if L<=R]
    K=[[int(one_cycle(caps[i],caps[j])) for j in right_ids] for i in left_ids]
    I,J=full_rank_minor(K)
    up_left=tuple(map(reverse,tails(a,t)))
    low_left=tuple(map(reverse,tails(b,t)))
    up_right=tails(a,R)
    low_right=tails(b,R)
    selected_left=[left_ids[i] for i in I]
    selected_right=[right_ids[j] for j in J]
    measures=lambda ids:[make_cap(a,b,caps[i]) for i in ids]
    left=closed_rows(up_left,low_left,measures(selected_right),False)
    right=closed_rows(up_right,low_right,measures(selected_left),True)
    inv=inverse([[K[i][j] for j in J] for i in I])
    for row in K:
        c=[sum(row[j]*v for j,v in zip(J,col)) for col in zip(*inv)]
        assert all(sum(x*K[i][q] for x,i in zip(c,I))==row[q] for q in range(len(row)))
    maximum_aux=max([(t+lengths[i])//2 for i in selected_right]+[(R+lengths[i])//2 for i in selected_left],default=0)
    assert maximum_aux<=(t+R)//2
    restricted=combine(left,right,inv,False)
    full_check=(a+b<=8)
    if full_check:
        all_caps=measures(range(len(caps)))
        lf=closed_rows(up_left,low_left,all_caps,False)
        rg=closed_rows(up_right,low_right,all_caps,True)
        full=combine(lf,rg,inverse_gram(a+b),True)
        assert full==restricted
    plain=sum(restricted.values())
    weighted=sum(centroid(p+u)*v for (p,u),v in restricted.items())
    return [a,b,len(caps),len(left_ids),len(right_ids),len(I),maximum_aux,plain,weighted,full_check]


def test_order(n,t):
    R=2*n-t
    literal,weighted_literal,checks=literal_totals(n,t)
    counts=[]
    for a in range(t%2,t+1,2):
        for b in range(t%2,t+1,2):
            row=by_height(a,b,t,R)
            assert tuple(row[7:9])==(literal[a,b],weighted_literal[a,b]),(n,t,a,b)
            counts.append(row)
    total=sum(r[7] for r in counts)
    assert total==sum(r[8] for r in counts)==sum(literal.values())==sum(weighted_literal.values())
    return dict(n=n,t=t,R=R,literal_word_pairs=checks,by_height=counts,closed_meander_number=total)


if __name__=="__main__":
    cli=ArgumentParser(description=__doc__)
    cli.add_argument("--write",type=Path)
    args=cli.parse_args()
    result=dict(source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),length_test=test_lengths(),
                complete_counts=[test_order(n,t) for n in range(1,9) for t in sorted({min(n,4),n})],
                columns=["top_height","lower_height","cap_count","left_support","right_support","compatible_rank",
                         "maximum_restricted_auxiliary_order","ordinary_contribution","centroid_contribution",
                         "full_inverse_comparison"])
    if args.write:
        pending=args.write.with_suffix(".pending")
        pending.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
        pending.replace(args.write)
    else:
        script=Path(__file__).resolve()
        default=script.parent/"certificates/all-boundaries.json"
        if not default.is_file():
            default=script.parents[3]/"output/research73/star/all-boundaries.json"
        assert result==json.loads(default.read_bytes())
    print(result["length_test"],[r["closed_meander_number"] for r in result["complete_counts"]])
