"""Exact Walsh joins on four transverse binary four-spaces, including the nonsemisimple case."""
from collections import Counter,defaultdict
import json
import numpy as np

from local import mul


def classify(L):
    x=np.arange(16)
    parity=np.array([i.bit_count()%2 for i in range(16)])
    H=1-2*parity[np.arange(1,16)[:,None]&x]
    signs=1-2*((np.arange(1<<15)[:,None]>>np.arange(15))&1)
    ts=signs@H
    ts=ts[(ts.min(axis=1)>=-6)&(ts.max(axis=1)<=10)]
    lookup={tuple(t%16):t for t in ts}
    rows=[t for t in ts if tuple(t)==min(tuple(t[x^s]) for s in x)]
    assert len(rows)*16==len(lookup)==7616
    u,v=np.meshgrid(x,x,indexing='ij')
    inv=np.argsort(x^L)
    assert np.array_equal(np.sort(x^L),x)
    Y=inv[u^v]; X=u^Y
    def defect(m):return ((m-m[0,:]-m[:,0,None]+m[0,0])%16).tobytes()
    right=defaultdict(list)
    for h in rows:right[defect(h[Y])].append(h)
    answers=[]
    for g in rows:
        for h in right[defect(-g[X])]:
            q=10-g[X]-h[Y]
            for t in range(16):
                c=lookup.get(tuple((q[:,0]-t)%16))
                b=lookup.get(tuple((q[0,:]-q[0,0]+t)%16))
                if c is None or b is None:continue
                assert set((g[X]+h[Y]+c[:,None]+b).flat)=={-6,10}
                profiles=[]
                for row in (g,h,c,b):
                    neg=(row@H.T)//16
                    assert set(neg.tolist())<={-1,1}
                    profiles.append(sum(1<<i for i,val in enumerate(neg) if val==-1))
                answers.append(profiles)
    signatures=Counter(tuple(sorted(p.bit_count()%2 for p in ps)) for ps in answers)
    return dict(translated_classes=len(answers),sign_parities={str(k):v for k,v in signatures.items()},
                representatives=answers)


if __name__=='__main__':
    # Companions of x^2+x+1 on two bit pairs, with a single Jordan coupling.
    def omega(a):return (0,2,3,1)[a]
    maps={'field_order3':np.array([mul(x,6) for x in range(16)]),
          'field_primitive':np.array([mul(x,2) for x in range(16)]),
          'jordan6':np.array([(omega(x%4)^(x//4))+4*omega(x//4) for x in range(16)])}
    print(json.dumps({lab:classify(L) for lab,L in maps.items()},indent=2))
