"""Construct a three-fibre Cayley control with a fixed-point-free order-three action."""
import argparse
from collections import Counter
from itertools import combinations,product
import json
from pathlib import Path
import sys

import numpy as np

from assemble import Ports,HERE
from construct768 import check,graph6,local_four_cliques


def search(planes):
    for name,lines in planes.items():
        if name=='field':continue
        P=Ports(lines)
        labels=np.full(256,-1)
        for j,line in enumerate(lines):labels[line]=j
        for i,j,k in combinations(range(17),3):
            a,b,c=(P.lines[t] for t in (i,j,k))
            coords={x^y:(ix,iy) for ix,x in enumerate(a) for iy,y in enumerate(b)}
            convert=dict(coords[u] for u in c)
            b=[b[convert[t]] for t in range(16)]
            # In these coordinates the three selected lines are (x,0),(0,y),(z,z).
            M=np.empty(256,dtype=int)
            for x,y in product(range(16),repeat=2):M[a[x]^b[y]]=a[y]^b[x^y]
            targets=labels[M[np.array(lines)]]
            if np.any(targets!=targets[:,0,None]):continue
            g=targets[:,0]
            cycles=[(l,g[l],g[g[l]]) for l in range(17) if l<min(g[l],g[g[l]])]
            for groups in combinations(cycles,4):
                for sel in product(*groups):
                    port=tuple(sorted(map(int,sel)))
                    val=P.edge(port)
                    if not val or val[0]!='jordan6':continue
                    raw=val[1][0][1]
                    z={p:(-1 if mask>>(v-1)&1 else 1) for idx,mask in zip(port,raw)
                       for v,p in enumerate(P.lines[idx]) if v}
                    assert np.array_equal(M[M[M]],np.arange(256)) and np.count_nonzero(M==np.arange(256))==1
                    return name,port,z,M
        P.edge.cache_clear()
    raise ValueError('No orbit of three disjoint Jordan supports was found')


def construct(planes):
    name,port,z,M=search(planes)
    zs=[z,{int(M[p]):s for p,s in z.items()},{int(M[M[p]]):s for p,s in z.items()}]
    assert len(set.union(*(set(v) for v in zs)))==180
    H=np.array([[1-2*((x&y).bit_count()%2) for x in range(256)] for y in range(256)])
    within=[];across=[]
    for i,edge in enumerate(zs):
        values=H[:,list(edge)]@np.array(list(edge.values()))
        assert Counter(values)=={-6:160,10:96}
        across.append(np.flatnonzero(values==10).tolist())
        active=list(set(edge)|set(zs[i-1]))
        vals=H[:,active].sum(axis=1)
        assert Counter(vals)=={120:1,8:120,-8:135}
        within.append(np.flatnonzero(vals==8).tolist())
    rows=[sum(1<<(256*j+(x^d)) for j,ds in ((i,within[i]),((i+1)%3,across[i]),((i-1)%3,across[i-1])) for d in ds)
          for i in range(3) for x in range(256)]
    from contextlib import redirect_stdout
    with redirect_stdout(sys.stderr):check(rows,(768,312,136,120))
    # The inverse transpose on point labels is dual to the frequency action.
    R=[sum(((int(M[M[1<<j]])&x).bit_count()%2)<<j for j in range(8)) for x in range(256)]
    maps=[list(range(256)),R,[R[R[x]] for x in range(256)]]
    D=[divmod(v,256) for v in range(768) if rows[0]>>v&1]
    for i,x in product(range(3),range(256)):
        assert rows[256*i+x]==sum(1<<(256*((i+j)%3)+(x^maps[i][y])) for j,y in D)
    result=dict(plane=name,seed_lines=port,negative_frequencies=sorted(p for p,s in z.items() if s<0),
                positive_frequencies=sorted(p for p,s in z.items() if s>0),
                frequency_basis_images=M[[1<<j for j in range(8)]].tolist(),
                vertex_basis_images=[R[1<<j] for j in range(8)],
                fixed_points=sum(i==r for i,r in enumerate(R)),
                four_cliques_per_fibre=[local_four_cliques(rows,256*i) for i in range(3)])
    return result,rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('graph',type=Path)
    args=p.parse_args()
    result,rows=construct(json.loads((HERE/'planes.json').read_text()))
    args.graph.write_text(graph6(rows))
    print(json.dumps(result,indent=2))
