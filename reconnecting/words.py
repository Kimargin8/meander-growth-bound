from functools import cache
from kernel import literal
@cache
def matchings(k):
    if not k:return [()]
    out=[]
    for j in range(k):
        b=1+2*j
        for inside in matchings(j):
            for tail in matchings(k-1-j):
                out.append((b,*(1+i for i in inside),0,*(b+1+i for i in tail)))
    return out


def open_pair(p,u):
    return tuple([i+(2 if i>=u else 0) for i in p[:u]]+[u+1,u]+
                 [i+(2 if i>=u else 0) for i in p[u:]])


def splice(p,u):
    if u==0 or u==len(p) or p[u-1]==u:return None
    adj=list(p)
    a,b=p[u-1],p[u]
    adj[a],adj[b]=b,a
    return tuple(i-(2 if i>u else 0) for i in adj[:u-1]+adj[u+1:])


def word(a,b):
    return "".join({(1,1):"U",(-1,-1):"D",(1,-1):"R",(-1,1):"L"}[
        (1 if a[i]>i else -1,1 if b[i]>i else -1)] for i in range(len(a)))


def accept(text,w,h):
    start=[];end={};height=[0];gap=[0]
    for i,c in enumerate(text):
        if c=="U":start.append(i)
        if c=="D":end[start.pop()]=i+1
        height.append(height[-1]+2*(int(c=="U")-int(c=="D")))
        gap.append(gap[-1]+(1 if c in "UR" else -1))
    def unit_children(i,j):
        for k in range(i+1,j-1):
            if text[k]=="U" and abs(gap[end[k]]-gap[k])!=1:return False
        return True
    def fixed(i,j,d):
        opens=[k for k in range(i+1,j-1) if text[k]=="U"]
        if len(opens)!=1:return False
        k,=opens;e=end[k]
        return (gap[k]==gap[i]+1 and gap[e]==gap[k]+d and
                max(k-i-1,e-k-2,j-e-1)<=6)
    pos=1;pending=False
    while pos<len(text)-1:
        if height[pos]!=w or text[pos]!="U":
            pending=False;pos+=1;continue
        j=end[pos];d=gap[j]-gap[pos]
        regular=abs(d) in (1,2,3) and max(height[pos:j+1])<=w+2*h and unit_children(pos,j)
        special=abs(d)==2 and fixed(pos,j,d)
        if not (regular or special) or (regular and special):return False
        if pending and d%2==0:return False
        pending=d%2==0
        pos=j
    return not pending


def direct(w,h,n):
    total=0
    for i,a in enumerate(matchings(n)):
        for j,b in enumerate(matchings(n)[:i+1]):
            p=b[a[0]];length=1
            while p:
                p=b[a[p]];length+=1
            if length==n and accept(word(a,b),w,h): total+=1+(i!=j)
    return total


def simple_walks(w,a,b):
    p=[0]*(w+1);p[a]=1;v=[]
    for n in range(7):
        v.append(p[b])
        p=[(p[i-1] if i else 0)+(p[i+1] if i<w else 0) for i in range(w+1)]
    return v


def computed(w,h,through):
    degree=2*through-2
    weights=literal(w,h,degree)
    fixed=[]
    for u in range(w-1):
        coeff=[0]*(degree+1)
        p=simple_walks(w+2,u+1,u+1)
        q=simple_walks(w+4,u+2,u+4)
        r=simple_walks(w+2,u+3,u+3)
        for a,A in enumerate(p):
            for b,B in enumerate(q):
                for c,C in enumerate(r):
                    if 4+a+b+c<=degree:coeff[4+a+b+c]+=A*B*C
        fixed.append(coeff)
    states=[(p,u) for k in range(1,w//2+1) for p in matchings(k) for u in range(2*k+1)]
    idx={s:i for i,s in enumerate(states)}
    odd=[[] for s in states];even=[[] for s in states]
    for i,(p,u) in enumerate(states):
        for g in (u-1,u+1):
            if 0<=g<=len(p):odd[i].append((1,idx[p,g],1))
        closed=splice(p,u)
        if closed is not None:odd[i].append((1,idx[closed,u-1],1))
        opened=open_pair(p,u)
        if len(p)<w:
            odd[i].append((1,idx[opened,u+1],1)); continue
        for delta in (-3,-2,-1,1,2,3):
            g=u+delta
            if not 0<=g<=w:continue
            d=abs(delta);small=min(u,g)
            if d==2:even[i].extend((k,idx[p,g],v) for k,v in enumerate(fixed[small]) if v)
            close=splice(opened,g+1)
            if close is not None:
                (even if d==2 else odd)[i].extend((k,idx[close,g],v) for k,v in enumerate(weights[d][small]) if v)
    init=[int(s==((1,0),1)) for s in states]
    normal=[init];waiting=[[0]*len(states)]
    for n in range(1,degree+1):
        o=[sum(v*normal[n-k][j] for k,j,v in row if k<=n) for row in odd]
        full=[o[i]+sum(v*waiting[n-k][j] for k,j,v in row if k<=n) for i,row in enumerate(even)]
        normal.append(full);waiting.append(o)
    return [normal[2*n-2][idx[((1,0),1)]] for n in range(1,through+1)]


