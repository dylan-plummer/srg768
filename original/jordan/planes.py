"""Generate order-16 spreads by binary semifield multiplication and net replacement."""
from collections import Counter,defaultdict
from itertools import combinations
import json
import sys

import numpy as np

from local import spaces,edge_type
from assemble import spreads,gl4


def fingerprint(lines):
    # For each four-set record whether its cross-ratio is semisimple cubic,
    # nonsemisimple, or irreducible. Aggregate at each line to distinguish classes.
    masks=[sum(1<<p for p in l) for l in lines]
    pointwise=[Counter() for _ in lines]
    for four in combinations(range(17),4):
        typ=edge_type([masks[i] for i in four])
        for i in four:pointwise[i][typ]+=1
    pivots={}
    def insert(row):
        while row:
            j=row.bit_length()-1
            if j not in pivots:
                pivots[j]=row
                return
            row^=pivots[j]
    insert(sum(1<<(256+i) for i in range(17)))
    for i,line in enumerate(lines):
        done=set()
        for x in range(256):
            if x in done:continue
            coset=[x]+[x^v for v in line]
            done.update(coset)
            insert((1<<(256+i))+sum(1<<v for v in coset))
    return len(pivots),tuple(sorted(tuple(c[t] for t in ('order3','jordan6')) for c in pointwise))


def semifields():
    gl=list(gl4())
    encode=lambda g:sum(g[1<<j]<<(4*j) for j in range(4))
    lookup={encode(g):g for g in gl}
    valid=np.zeros(65536,dtype=bool);valid[list(lookup)]=True
    ids=np.array(sorted(lookup))
    eye=encode(list(range(16)))
    rng=np.random.default_rng(1024)
    for a in rng.permutation(ids[valid[ids^eye]]):
        left=[0,eye,int(a),eye^int(a)]
        B=ids[np.all(valid[ids[:,None]^left],axis=1)]
        for b in rng.choice(B,size=min(8,len(B)),replace=False):
            base=left+[v^int(b) for v in left]
            C=ids[np.all(valid[ids[:,None]^base],axis=1)]
            if len(C):
                c=int(rng.choice(C))
                mats=base+[v^c for v in base]
                tables=[[0]*16 if not k else lookup[k] for k in mats]
                yield [[x+16*g[x] for x in range(1,16)] for g in tables]+[[16*y for y in range(1,16)]]


def replacements(lines,subspace_points):
    labels=np.zeros(256,dtype=np.int32)
    for i,l in enumerate(lines):labels[l]=1<<i
    footprint=np.bitwise_or.reduce(labels[subspace_points],axis=1)
    masks=[sum(1<<p for p in l) for l in lines]
    small=defaultdict(list)
    for i in np.flatnonzero(np.isin(footprint,[x for x in np.unique(footprint) if int(x).bit_count()==5])):
        small[int(footprint[i])].append(subspace_points[i].tolist())
    for occupied,family in small.items():
        if len(family)!=5:continue
        if len(set(sum(family,[])))==75:
            yield [l for i,l in enumerate(lines) if not (occupied>>i)&1]+family


def main():
    points=np.array(list(spaces()),dtype=np.uint8)
    classes={}
    def add(name,lines):
        sig=fingerprint(lines)
        if sig in classes:return
        assert len(set(sum(lines,[])))==255
        classes[sig]=(name,lines)
        print(name,sig[0],Counter(sig[1]),file=sys.stderr,flush=True)
    for name,lines in spreads():add(name,lines)
    # Eight is the published number of translation planes, although fingerprints
    # need not distinguish every pair. Keep the construction search bounded.
    for i,lines in enumerate(semifields()):
        add('semifield_'+str(i),lines)
        if len(classes)>=4 or i>=150:break
    seeds=list(classes.values())
    for name,lines in seeds:
        for j,new in enumerate(replacements(lines,points)):add(name+'_derived_'+str(j),new)
    for name,lines in list(classes.values()):
        dual=[[v for v in range(1,256) if all((v&x).bit_count()%2==0 for x in line)] for line in lines]
        add(name+'_perp',dual)
    print(json.dumps({name:lines for name,lines in classes.values()}))


if __name__=='__main__':main()
