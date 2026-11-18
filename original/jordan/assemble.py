"""Assemble signed 256-point cross blocks on a binary partial spread."""
from collections import Counter
from functools import lru_cache
from itertools import combinations,product
import json
from pathlib import Path
import sys

import numpy as np

from local import bits,mul,edge_type

HERE=Path(__file__).resolve().parent
CLASSES=json.loads((HERE/'crossratio.json').read_text())
PAIRS=tuple(combinations(range(4),2))


def span(basis):
    vals=[0]
    for b in basis:vals += [v^b for v in vals]
    return vals


def coordinate_space(line):
    basis=[]
    s={0}
    for v in sorted(line):
        if v not in s:
            basis.append(v);s.update([b^v for b in s])
    assert len(basis)==4
    return span(basis)


def spreads():
    field=[[x+16*mul(m,x) for x in range(1,16)] for m in range(16)]+[[16*x for x in range(1,16)]]
    yield 'field',field
    K=[0,1,6,7]
    def opposite(a):
        return [mul(a,x)+16*mul(a,y) for x,y in product(K,repeat=2) if x or y]
    lines=[field[m] for m in range(16) if m not in K]
    done=set()
    for a in range(1,16):
        if a not in done:
            done.update(mul(a,i) for i in K)
            lines.append(opposite(a))
    assert len(lines)==17 and len(set(sum(lines,[])))==255
    yield 'hall',lines


def gl4():
    for a,b in combinations(range(1,16),2):
        for x,y in ((a,b),(b,a)):
            for c in set(range(1,16))-{x,y,x^y}:
                for d in set(range(1,16))-set(span([x,y,c])):
                    yield span([x,y,c,d])


@lru_cache(None)
def shared_tables():
    gl=np.array(list(gl4()),dtype=np.uint8)
    chars=[sum(((k&i).bit_count()%2)<<(i-1) for i in range(1,16)) for k in range(16)]
    return gl,chars,[min(i^t for t in chars) for i in range(1<<15)]


@lru_cache(None)
def intertwiner(target,source):
    gl=shared_tables()[0]
    valid=np.flatnonzero(np.all(np.array(target)[gl]==gl[:,source],axis=1))
    return None if not len(valid) else gl[valid[0]]


class Ports:
    def __init__(self,lines):
        self.lines=[coordinate_space(l) for l in lines]
        self.gl,self.chars,self.reduce=shared_tables()
        L0=[mul(x,6) for x in range(16)]
        w=[0,2,3,1]
        L1=[(w[x%4]^(x//4))+4*w[x//4] for x in range(16)]
        self.types={name:[sum(((m&L[1<<i]).bit_count()%2)<<i for i in range(4))
                          for m in range(16)] for name,L in [('field_order3',L0),('jordan6',L1)]}

    @lru_cache(None)
    def edge(self,four):
        a,b,c,d=[self.lines[i] for i in four]
        coords={x^y:(i,j) for i,x in enumerate(a) for j,y in enumerate(b)}
        pc=dict(coords[x] for x in c)
        pd=dict(coords[x] for x in d)
        cb={v:k for k,v in pc.items()}
        L=np.array([cb[pd[i]] for i in range(16)])
        # The fourth cross ratio must be one of the two types with a Boolean signing.
        for typ,R in self.types.items():
            if typ=='field_order3' and not np.array_equal(L[L]^L,np.arange(16)):continue
            if typ=='jordan6' and not np.array_equal(L[L[L[L[L[L]]]]],np.arange(16)):continue
            f=intertwiner(tuple(L),tuple(R))
            if f is None:continue
            phi=[(a[f[i]],b[pc[f[i]]],a[f[i]]^b[pc[f[i]]],a[f[i]]^b[pc[f[R[i]]]]) for i in range(16)]
            perm=[[self.lines[k].index(z[j])-1 for z in phi[1:]] for j,k in enumerate(four)]
            opts=[]
            for signed in CLASSES[typ]['representatives']:
                vals=tuple(sum(1<<perm[j][i] for i in bits(mask)) for j,mask in enumerate(signed))
                opts.append((tuple(self.reduce[x] for x in vals),vals))
            return typ,opts
        return None

    def pair(self,a,b):
        four=tuple(sorted(a+b))
        v=self.edge(four)
        if v is None:return None
        indices=[four.index(x) for x in a+b]
        return [(tuple(red[i] for i in indices),tuple(raw[i] for i in indices)) for red,raw in v[1]]

    def solve(self):
        pairs=list(combinations(range(len(self.lines)),2))
        viable={a:set() for a in pairs}
        for a,b in combinations(pairs,2):
            if not set(a)&set(b) and self.edge(tuple(sorted(a+b))):
                viable[a].add(b);viable[b].add(a)
        stats=Counter()
        @lru_cache(None)
        def opts(a,b):return self.pair(a,b)
        @lru_cache(None)
        def can_meet(a,b,c,d):
            first={x[0][:2] for x in opts(a,b)}
            second={x[0][:2] for x in opts(a,c)}
            third={x[0][:2] for x in opts(a,d)}
            return any((p[0]^q[0],p[1]^q[1]) in third for p in first for q in second)
        # Canonical reduction modulo addition of characters is linear.
        assert all(self.reduce[v^t]==self.reduce[v] for v in range(1<<15) for t in self.chars)
        for a in pairs:
            for b in sorted(x for x in viable[a] if x>a):
                shared=viable[a]&viable[b]
                for c in sorted(x for x in shared if x>b):
                    for d in (x for x in shared&viable[c] if x>c):
                        layout=(a,b,c,d)
                        stats['layouts']+=1
                        if not all(can_meet(p,*(q for q in layout if q!=p)) for p in layout):continue
                        stats['four_local_checks']+=1
                        hit=self.join(layout,opts,stats)
                        if hit:return stats,hit
        return stats,None

    def join(self,layout,opts,stats):
        # With one line pair per omitted vertex, an edge is indexed by the two faces that contain it.
        choices=[opts(layout[i],layout[j]) for i,j in PAIRS]
        by=[{} for _ in choices]
        for e,options in enumerate(choices):
            for reduced,raw in options:by[e].setdefault(reduced[:2],[]).append((reduced,raw))
        def dfs(k,acc,used):
            if k==6:
                stats['signed_mod_linear']+=1
                return self.translations(layout,used)
            i,j=PAIRS[k]
            source=by[k].get(acc[i],[]) if (i,j) in [(0,3),(1,3),(2,3)] else choices[k]
            for reduced,raw in source:
                v=list(acc)
                v[i]=tuple(x^y for x,y in zip(v[i],reduced[:2]))
                v[j]=tuple(x^y for x,y in zip(v[j],reduced[2:]))
                if k==5 and (v[i]!=(0,0) or v[j]!=(0,0)):continue
                hit=dfs(k+1,v,used+[raw])
                if hit:return hit
            return None
        return dfs(0,[(0,0)]*4,[])

    def translations(self,layout,used):
        sums=[[0,0] for _ in layout]
        for (i,j),v in zip(PAIRS,used):
            for u,r in [(i,v[:2]),(j,v[2:])]:
                for k in range(2):sums[u][k]^=r[k]
        defects=[]
        for pair,vals in zip(layout,sums):
            opts=[t for t in range(256) if all(sum(((x&t).bit_count()%2)<<(h-1)
                  for h,x in enumerate(self.lines[l]) if h)==v for l,v in zip(pair,vals))]
            assert len(opts)==1
            defects+=opts
        if defects[0]^defects[1]^defects[2]^defects[3]:return None
        # Assign translations on a star to cancel all four defects.
        trans=[defects[j] if i==0 else 0 for i,j in PAIRS]
        return dict(layout=layout,signs=used,translations=trans,lines=self.lines)


if __name__=='__main__':
    out={}
    for name,lines in spreads():
        search=Ports(lines);stats,hit=search.solve()
        out[name]=dict(counts=dict(stats),hit=hit)
        print(name,dict(stats),bool(hit),file=sys.stderr,flush=True)
    print(json.dumps(out,indent=2))
