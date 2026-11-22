"""Граница восстановления взвешенной склейки по более коротким замыканиям.

Для заданного сектора (a,b) упорядочим верхние концы от старшего к младшему,
нижние — обратно. K[p,q]=1, когда наложение шапок даёт один цикл. Половины
имеют единичные веса; x[p] — число допустимых пар префиксов настоящей длины n.
I состоит из шапок минимальной длины L=a+b-j<n, E из шапок с L=n; j есть
число верхне-нижних соединений. Измерение каждой дешёвой шапкой есть (Kx)_I.
При известных x_E эти измерения определяют A x_I, где A=K_II. Слагаемое
x_I^T A x_I постоянно на каждом слое этой линейной системы. Полная симметричная
склейка x^T K x постоянна на всех слоях в точности тогда, когда
f=K_IE x_E ортогонален ker A: замена x_I на x_I+z, Az=0, меняет сумму на
2 f^T z. Для двух разных сторон условие ставится для каждой требуемой
линейной функционали отдельно. Это критерий восстановления из указанных
измерений: значения реальных масс могут подчиняться дополнительным тождествам.

Бесконечная помеха. При n=2k>=4 и положительных чётных a+b=2k шапка
A=(01)(23)... имеет j=0, а циклически сдвинутая B=(12)(34)...(2k-1,0)
имеет j=2. В полном K единственный партнёр B — A. Действительно, стянем
соседние пары B к k точкам на окружности. Связная вторая шапка даёт
некрещивающийся цикл через все k точек. Для k>=3 он идёт в граничном порядке:
любая диагональ разделила бы остающийся путь на две стороны. У каждой
соседней пары B только одна расстановка двух концов сохраняет отсутствие
пересечений с предыдущей и следующей парами. При k=2 проверка непосредственная.
Так B даёт нулевую строку K_II, а f_B=x_A>0. В минимальной половине A каждая
пара верхних либо нижних концов занимает две соседние вершины: между ними
не может находиться конец противоположного банка, а концов собственного
банка между соседями A нет. Перемешивание этих k блоков даёт x_A=binom(k,a/2).

В секторе (a,b)=(2k-2,2) указанное препятствие имеет явное решение x_B=k^2.
Единственный верхний закрывающий шаг образует дугу между соседними вершинами.
Удалим их, соединив два различных нижних ребра, к которым они принадлежали.
У полученной шапки B k дуг между циклическими соседями. Для обратной вставки
выбираем один из 2k-1 промежутков на верхнем банке и дугу, инцидентную грани
этого промежутка: разделяем дугу и присоединяем её концы к двум новым вершинам.
Некрещивающееся соединение единственно; дуга вне этой грани недоступна.
Из промежутков k лежат в однодуговых гранях, а k-1 лежат в центральной грани,
инцидентной всем k дугам. Получаем k+(k-1)k=k^2. Вклад A и B в сектор равен
2k^3. При n=4 настоящие внутренние массы (4) и крайние (2), вклад равен 16.
При n=8 индивидуальное знание всех изолированных строк само по себе уже не
чинит ранговое восстановление. При n=10 этого не дают совместно даже маргиналы
по (L, нараяновскому рангу) и знание изолированных строк.

Также проверяется конечный положительный вариант: если задать каждую массу
с L=n-2, у усиленной матрицы наблюдений во всех проверенных смешанных секторах
n<=12 полный столбцовый ранг. Это конечная проверка; оценка числа либо цены
нужных замыканий и утверждение для всех n из неё не следуют.

Проверка ограниченного избытка. У сквозного пути минимальная длина 1, у
пути с концами на одной стороне — 2. Большая длина отличается на чётное
число. При суммарном избытке 2d удаление целых растянутых путей затрагивает
не более d компонент и 4d вершин. Оставшиеся вершины имеют хотя бы один
внешний конец каждая. Его номер в соответствующем банке однозначно задаёт
положение вершины относительно каждой марки, если даны числа более ранних
обычных концов в обоих банках. Поэтому пересечения между маркированными и
обычными путями, в том числе с лучами к краю, полностью определены метками.
Обычные сквозные вершины продвигают оба банка одновременно. Вставка события
нижнего банка между двумя концами открытой верхне-верхней пары запрещена:
она пересекла бы её соединительное нижнее ребро. Для другого банка действует
симметричное правило. Динамика на (i,j) имеет не более трёх переходов и
считает перестановки минимальной половины между последовательными марками.
Здесь независимо проверены все описания, имеющие хотя бы одно настоящее
продолжение, до n=10 при d=1,2. Число их истинных реализаций сравнивается
как с полным списком минимальных половин, так и с этой двухбанковой динамикой.
Для реализации общего алгоритма дополнительно требуется проверять допустимость
самих помеченных дуг, границу шапки и все пересечения с внешними лучами.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb, gcd, lcm
from hashlib import sha256
from pathlib import Path
import json
import sys


@lru_cache(None)
def half_words(length):
    """A side of a prefix, with unclosed arrows listed from oldest to newest."""
    states=[((),())]
    for i in range(length):
        nxt=[]
        for edges,stack in states:
            nxt.append((edges,stack+(i,)))
            if stack:
                nxt.append((edges+((stack[-1],i),),stack[:-1]))
        states=nxt
    return tuple(states)


@lru_cache(None)
def populations(n):
    by_shape=defaultdict(Counter)
    for top,upper in half_words(n):
        for bottom,lower in half_words(n):
            if not upper and not lower:
                continue
            parent=list(range(n))
            def find(i):
                while parent[i]!=i:
                    i=parent[i]
                return i
            for x,y in top+bottom:
                p,q=find(x),find(y)
                if p==q:
                    break  # a closed cycle has already separated
                parent[p]=q
            else:
                ports=upper+lower[::-1]
                paths=defaultdict(list)
                for i,v in enumerate(ports):
                    paths[find(v)].append(i)
                assert set(paths)=={find(i) for i in range(n)} and all(len(v)==2 for v in paths.values())
                mates=[-1]*len(ports)
                for i,j in paths.values():
                    mates[i],mates[j]=j,i
                by_shape[len(upper),len(lower)][tuple(mates)]+=1
    return by_shape


def connected(p,q):
    v,n=0,len(p)//2
    for i in range(1,n+1):
        v=p[q[v]]
        if v==0:
            return int(i==n)
    raise AssertionError("a cycle must return")


def nullspace(A):
    rows=[list(map(Fraction,r)) for r in A]
    pivots=[]
    for j in range(len(rows[0])):
        i=len(pivots)
        pivot=next((k for k in range(i,len(A)) if rows[k][j]),None)
        if pivot is None:
            continue
        rows[i],rows[pivot]=rows[pivot],rows[i]
        m=rows[i][j]
        rows[i]=[v/m for v in rows[i]]
        for k,r in enumerate(rows):
            if k!=i and r[j]:
                scale=r[j]
                rows[k]=[v-scale*w for v,w in zip(r,rows[i])]
        pivots.append(j)
    basis=[]
    for j in range(len(rows[0])):
        if j not in pivots:
            v=[Fraction(0)]*len(rows[0])
            v[j]=1
            for i,p in enumerate(pivots):
                v[p]=-rows[i][j]
            basis.append(v)
    return basis


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def grade(p):
    q=[p[2*i+1]//2 for i in range(len(p)//2)]
    groups=0
    remaining=set(range(len(q)))
    while remaining:
        i=next(iter(remaining))
        groups+=1
        while i in remaining:
            remaining.remove(i)
            i=q[i]
    return groups


def scan(limit,controls=()):
    singular=0
    examined=0
    for n in range(1,limit+1):
        data=populations(n)
        for (a,b),masses in sorted(data.items()):
            if a<b:
                continue
            def cost(p):
                return a+b-sum(s>=a for s in p[:a])
            interior=sorted(p for p in masses if cost(p)<n)
            exterior=sorted(p for p in masses if cost(p)==n)
            assert len(interior)+len(exterior)==len(masses)
            if not interior or not exterior:
                continue
            examined+=1
            A=[[connected(p,q) for q in interior] for p in interior]
            observations=A.copy()
            if "total" in controls:
                observations.append([1]*len(interior))
            if "isolated" in controls:
                for i,row in enumerate(A):
                    if not any(row):
                        observations.append([int(k==i) for k in range(len(interior))])
            if "near" in controls:
                for i,p in enumerate(interior):
                    if cost(p)==n-2:
                        observations.append([int(k==i) for k in range(len(interior))])
            if "grades" in controls:
                stats=[(cost(p),grade(p)) for p in interior]
                observations.extend([int(x==y) for x in stats] for y in sorted(set(stats)))
            if "symmetry" in controls and a==b:
                for i,p in enumerate(interior):
                    q=tuple(len(p)-1-v for v in p[::-1])
                    j=interior.index(q)
                    assert masses[p]==masses[q]
                    if i<j:
                        observations.append([int(k==i)-int(k==j) for k in range(len(interior))])
            kernel=nullspace(observations)
            if not kernel:
                continue
            singular+=1
            ext=[sum(connected(p,q)*masses[q] for q in exterior) for p in interior]
            for z in kernel:
                effect=dot(z,ext)
                if effect:
                    den=lcm(*(v.denominator for v in z))
                    ints=[int(v*den) for v in z]
                    div=gcd(*ints)
                    z=[v//div for v in ints]
                    effect=dot(z,ext)
                    assert all(dot(row,z)==0 for row in observations)
                    clipped1=[masses[p]-max(0,-v) for p,v in zip(interior,z)]
                    clipped2=[masses[p]-max(0,v) for p,v in zip(interior,z)]
                    assert min(clipped1+clipped2)>=0
                    x=[masses[p] for p in interior]
                    E=[[connected(p,q) for q in exterior] for p in exterior]
                    y=[masses[p] for p in exterior]
                    tail=dot(y,[dot(r,y) for r in E])
                    value=lambda xx:dot(xx,[dot(r,xx) for r in A])+2*dot(xx,ext)+tail
                    assert value(clipped1)-value(clipped2)==2*effect
                    return dict(n=n,a=a,b=b,interior=interior,exterior=exterior,
                                mass_interior=x,mass_exterior=y,controls=controls,full_sector=value(x),
                                gram_interior=A,exterior_functional=ext,
                                kernel_direction=z,effect=effect,actual_subfamilies=[clipped1,clipped2],
                                subfamily_sectors=[value(clipped1),value(clipped2)],singular_checked=singular,
                                mixed_sector_representatives_checked=examined,gram_rank=len(A)-len(nullspace(A)))
    return dict(max_n=limit,singular_checked=singular,mixed_sector_representatives_checked=examined)


def marked_description(n,top,bottom):
    """Trace complete components; every retained ordinary vertex has a named external end."""
    upper_edges,upper_ends=top
    lower_edges,lower_ends=bottom
    a,b=len(upper_ends),len(lower_ends)
    edges=[[-1]*n,[-1]*n]
    for side,group in enumerate((upper_edges,lower_edges)):
        for i,j in group:
            edges[side][i],edges[side][j]=j,i
    for port,i in enumerate(upper_ends+lower_ends[::-1]):
        edges[int(port>=a)][i]=n+port
    cap=[-1]*(a+b)
    marked=set()
    covered=set()
    stretches=0
    for p,start in enumerate(upper_ends+lower_ends[::-1]):
        if cap[p]>=0:
            continue
        side=1-int(p>=a)
        i=start;vertices=[i]
        while edges[side][i]<n:
            i=edges[side][i]
            vertices.append(i)
            side=1-side
        end=edges[side][i]-n
        cap[p],cap[end]=end,p
        assert not covered.intersection(vertices)
        covered.update(vertices)
        shortest=2 if (p<a)==(end<a) else 1
        excess=len(vertices)-shortest
        assert excess>=0 and excess%2==0
        if excess:
            stretches+=1
            marked.update(vertices)
    if len(covered)!=n:
        return None
    L=a+b-sum(end>=a for end in cap[:a])
    d=(n-L)//2
    assert stretches<=d and len(marked)<=4*d
    ordinary_ports=tuple(p for p,v in enumerate(upper_ends+lower_ends[::-1]) if v not in marked)
    base=(sum(p<a for p in ordinary_ports),sum(p>=a for p in ordinary_ports),
          tuple(ordinary_ports.index(cap[p]) for p in ordinary_ports))
    new_id={p:i for i,p in enumerate(sorted(marked))}
    marks=[];progress=[];gaps=[(0,0)]
    for p in range(n):
        if p in marked:
            ends=tuple(new_id[q] if q<n else len(marked)+q-n for q in (edges[0][p],edges[1][p]))
            marks.append(ends)
            progress.append(gaps[-1])
        else:
            step=int(edges[0][p]>=n),int(edges[1][p]>=n)
            assert any(step)
            gaps.append((gaps[-1][0]+step[0],gaps[-1][1]+step[1]))
    desc=(a,b,tuple(cap),tuple(marks),tuple(progress))
    return d,base,desc,frozenset(gaps)


@lru_cache(None)
def number_minimal_paths(base,start,end):
    """Ordinary shuffles in an interval of the two ordered banks; through events advance both."""
    a,b,p=base
    p_lower=tuple(a+b-1-p[a+b-1-j] if p[a+b-1-j]>=a else None for j in range(b))
    def active_top(i):
        return any(s<i<=p[s]<a for s in range(i))
    def active_lower(j):
        return any(p_lower[s] is not None and j<=p_lower[s] for s in range(j))
    ways=Counter({start:1})
    for i in range(start[0],end[0]+1):
        for j in range(start[1],end[1]+1):
            v=ways[i,j]
            if not v:
                continue
            if i<end[0] and p[i]<a and not active_lower(j):
                ways[i+1,j]+=v
            if j<end[1] and p_lower[j] is not None and not active_top(i):
                ways[i,j+1]+=v
            if i<end[0] and j<end[1] and p[i]==a+b-1-j and not active_top(i) and not active_lower(j):
                ways[i+1,j+1]+=v
    return ways[end]


def marker_check(max_n=10):
    """For every actually liftable description, count all ordinary compatible shuffles."""
    minimal=defaultdict(list)
    decorated=Counter()
    base_for={}
    for n in range(max_n+1):
        for top in half_words(n):
            for bottom in half_words(n):
                info=marked_description(n,top,bottom)
                if info is None:
                    continue
                d,base,desc,gaps=info
                if not d:
                    minimal[base].append(gaps)
                elif d<=2:
                    key=(n,d,desc)
                    decorated[key]+=1
                    if key in base_for:
                        assert base_for[key]==base
                    base_for[key]=base
    counts=Counter()
    repeated=Counter()
    maximum=0
    for base,gap_sets in minimal.items():
        assert number_minimal_paths(base,(0,0),base[:2])==len(gap_sets)
    for key,true_count in decorated.items():
        n,d,desc=key
        allowed=set(desc[-1])
        base=base_for[key]
        actual=sum(allowed<=gaps for gaps in minimal[base])
        assert actual==true_count,(key,true_count,actual)
        stops=((0,0),)+desc[-1]+(base[:2],)
        dp_count=1
        for u,v in zip(stops,stops[1:]):
            dp_count*=number_minimal_paths(base,u,v)
        assert dp_count==true_count
        counts[d]+=1
        repeated[d]+=actual>1
        maximum=max(maximum,actual)
    return dict(max_n=max_n,excess_pairs=dict(counts),descriptions_with_multiple_shuffles=dict(repeated),
                true_nonminimal_halves=sum(decorated.values()),maximum_equivalent_shuffles=maximum)


def verify():
    root=Path(__file__).resolve().parents[3]
    known_path=root/"references/oeis/2026-11-15/b005315.txt"
    known={int(s):int(v) for s,v in (l.split() for l in known_path.read_text().splitlines() if l.strip())}
    assembled={}
    for n in range(1,13):
        reconstructed=0
        for (a,b),row in populations(n).items():
            if a<b:
                continue
            total=sum(m*v*connected(p,q) for p,m in row.items() for q,v in row.items())
            reconstructed+=total*(1+(a!=b))
        assert reconstructed==known[n]
        assembled[n]=reconstructed
    even=[]
    for k in range(2,7):
        n=2*k
        outer=tuple(i^1 for i in range(n))
        inner=tuple(((((i-1)%n)^1)+1)%n for i in range(n))
        for a in range(2,n,2):
            b=n-a
            row=populations(n)[a,b]
            assert row[outer]==comb(k,a//2)
            assert [p for p in row if connected(inner,p)]==[outer]
            if b==2:
                assert row[inner]==k*k
                even.append(dict(n=n,a=a,b=b,outer_mass=row[outer],inner_mass=row[inner],
                                 isolated_pair_contribution=2*row[inner]*row[outer]))
    controls=[(),("isolated",),("total","symmetry"),("isolated","grades","symmetry")]
    witnesses=[scan(10,c) for c in controls]
    assert [v["n"] for v in witnesses]==[4,8,8,10]
    assert [v["effect"] for v in witnesses]==[2,6,2,10]
    near=scan(12,("near",))
    assert near["singular_checked"]==0 and near["max_n"]==12
    return dict(convention="Upper endpoints oldest first, lower endpoints newest first. All counts are unweighted.",
                closed_counts=assembled,infinite_even_examples=even,witnesses=witnesses,
                all_near_layer_known=near,marker_progress=marker_check(),
                inputs={p.relative_to(root).as_posix():sha256(p.read_bytes()).hexdigest()
                        for p in (Path(__file__).resolve(),known_path)})


if __name__=="__main__":
    if not __debug__:
        raise ValueError("Run with Python assertions enabled")
    dest=Path(__file__).resolve().parents[3]/"output/research73/wholecut/strict-order-kernel.json"
    result=json.loads(json.dumps(verify()))
    if "--write" in sys.argv:
        dest.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    else:
        assert result==json.loads(dest.read_bytes())
    print("witness orders",[r["n"] for r in result["witnesses"]])
    print(result["all_near_layer_known"])
    print(result["marker_progress"])
