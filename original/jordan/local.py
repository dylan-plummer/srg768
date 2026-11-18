"""Enumerate four-space edge supports inside six Desarguesian spread members."""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import json
import sys

OLD = Path(__file__).resolve().parents[1]   # srg768/original (construct768.py); was ../construction_routes_20261118
sys.path.insert(0,str(OLD))
from construct768 import mul


def bits(S):
    while S:
        b=S&-S
        yield b.bit_length()-1
        S^=b


def slope_action():
    inv=[0]+[next(y for y in range(1,16) if mul(x,y)==1) for x in range(1,16)]
    return [[x^1 for x in range(16)]+[16], [mul(2,x) for x in range(16)]+[16],
            [16]+inv[1:]+[0], [mul(x,x) for x in range(16)]+[16]]


def representatives():
    gens=slope_action()
    unseen=set(combinations(range(17),6))
    while unseen:
        s=min(unseen)
        orbit={s}
        stack=[s]
        while stack:
            t=stack.pop()
            for g in gens:
                u=tuple(sorted(g[i] for i in t))
                if u not in orbit:orbit.add(u);stack.append(u)
        unseen-=orbit
        yield s,len(orbit)


def spaces():
    """Unique reduced binary bases, with the lowest set bit as pivot."""
    for pivots in combinations(range(8),4):
        cols=[(i,c) for i,p in enumerate(pivots) for c in range(p+1,8) if c not in pivots]
        for choices in range(1<<len(cols)):
            basis=[1<<p for p in pivots]
            for k in bits(choices):
                i,c=cols[k];basis[i]|=1<<c
            pts=[0]
            for b in basis:pts += [x^b for x in pts]
            yield pts[1:]


def edge_type(groups):
    """Coordinate-free characteristic of the cross-ratio of four disjoint spaces."""
    a,b,c,d=map(lambda x:[0]+list(bits(x)),groups)
    coords={u^v:(u,v) for u in a for v in b}
    assert len(coords)==256
    tc=dict(coords[p] for p in c)
    td=dict(coords[p] for p in d)
    back={v:u for u,v in tc.items()}
    z={u:back[td[u]] for u in a}
    if all(z[z[u]]^z[u]==u for u in a):return 'order3'
    if all(z[z[z[z[z[z[u]]]]]]==u for u in a):return 'jordan6'
    return 'irreducible4'


def local_catalogues():
    lines=[sum(1<<(a+16*mul(m,a)) for a in range(1,16)) for m in range(16)]
    lines.append(sum(1<<(16*b) for b in range(1,16)))
    slopes={p:i for i,l in enumerate(lines) for p in bits(l)}
    by=defaultdict(list)
    count=0
    for pts in spaces():
        count+=1
        occupied=frozenset(slopes[p] for p in pts)
        if len(occupied)<=6:by[occupied].append(sum(1<<p for p in pts))
    assert count==200787
    for rep,orbit in representatives():
        root=sum(lines[i] for i in rep)
        subs=[s for ls,vs in by.items() if ls<=set(rep) for s in vs]
        families=defaultdict(list)
        types=Counter()
        for four in combinations(subs,4):
            support=four[0]|four[1]|four[2]|four[3]
            if support.bit_count()!=60:continue
            typ=edge_type(four)
            types[typ]+=1
            if typ in ('order3','jordan6'):families[support].append(four)
        charts=[]
        for s,t in combinations(sorted(families),2):
            u=s^t
            if (s|t)==root and t<u and u in families:charts.append((s,t,u))
        yield dict(slopes=rep,orbit=orbit,subspaces=len(subs),decomposition_types=dict(types),
                   supports=len(families),charts=len(charts),mixed_charts=sum(
                   any(any(0<(s&l).bit_count()<15 for l in lines) for s in chart) for chart in charts)),families,charts,root,lines


if __name__=='__main__':
    print(json.dumps([r for r,*_ in local_catalogues()],indent=2))
