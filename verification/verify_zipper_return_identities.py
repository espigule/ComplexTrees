#!/usr/bin/env python3
"""Independent exact check of the seven endpoint-return rows.
Uses an integer-coefficient bivariate polynomial implementation, with q and Q
as independent symbols exchanged by conjugation. Polynomial
identities in these two symbols imply the physical identities Q=conj(q).
Words use the manuscript alphabet 1,2 and compose leftmost map outermost.
"""
import json
from pathlib import Path
from numbers import Integral
class Poly:
    def __init__(self,p=0):
        self.d={k:v for k,v in (p if isinstance(p,dict) else {(0,0):p}).items() if v}
    def __add__(self,o):
        o=o if isinstance(o,Poly) else Poly(o); d=self.d.copy()
        for k,v in o.d.items(): d[k]=d.get(k,0)+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self): return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,o): return self+-1*o
    def __rsub__(self,o): return o+-self
    def __mul__(self,o):
        o=o if isinstance(o,Poly) else Poly(o);d={}
        for (a,b),x in self.d.items():
            for (c,e),y in o.d.items():d[a+c,b+e]=d.get((a+c,b+e),0)+x*y
        return Poly(d)
    __rmul__=__mul__
    def __eq__(self,o): return self.d==(o.d if isinstance(o,Poly) else Poly(o).d)
    def __repr__(self): return str(self.d)
class s:
    Integer=Poly
    expand=staticmethod(lambda z:z)
    factor=staticmethod(lambda z:z)
q,Q=Poly({(1,0):1}),Poly({(0,1):1})
b=1-q

def bar(z):
    return Poly({(b,a):v for (a,b),v in z.d.items()})

def compose(f,g):
    a,t,k=f; A,T,K=g
    return s.expand(a*(bar(A) if k else A)),s.expand(t+a*(bar(T) if k else T)),bool(k)^bool(K)

def word(fs,w):
    f=(s.Integer(1),s.Integer(0),False)
    for c in w: f=compose(f,fs[int(c)-1])
    return f

rows=[
('DD01',(0,0,0,1),'12','22','1','1',q),
('DD10',(0,0,1,0),'11','21','2','2',b),
('DD11',(0,0,1,1),'1','2','12','21',q*b),
('DO10',(0,1,1,0),'11','21','22','22',b*(1-Q)),
('DO11',(0,1,1,1),'1','2','1212','2121',q*Q*b*(1-Q)),
('OO01',(1,1,0,1),'12','22','11','11',q*Q),
('OO10',(1,1,1,0),'11','21','22','22',b*(1-Q)),
]
results=[]
for name,(k1,k2,e1,e2),U,V,r,t,eta in rows:
    fs=((-q if e1 else q,q if e1 else s.Integer(0),k1),(-b if e2 else b,s.Integer(1) if e2 else q,k2))
    H=(s.expand(eta),s.expand(q-eta*q),False)
    defects=[]
    for W,R in [(U,r),(V,t)]:
        lhs=word(fs,W+R); rhs=compose(H,word(fs,W))
        defects += [s.expand(lhs[0]-rhs[0]),s.expand(lhs[1]-rhs[1]),int(lhs[2]!=rhs[2])]
    ok=all(v==0 for v in defects)
    assert ok,(name,defects)
    # Establish the seed maps really contain q as an endpoint image.
    seed_endpoints=[]
    for W in (U,V):
        a,z,k=word(fs,W)
        good=[j for j in [0,1] if s.expand(z+a*j-q)==0]
        assert good,(name,W)
        seed_endpoints.append(good)
    results.append(dict(family=name,seed=[U,V],returns=[r,t],eta=str(s.factor(eta)),identities_verified=ok,q_preimage_endpoints=seed_endpoints))
print(json.dumps(results,indent=2))
if __name__=='__main__':
    Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2)+'\n')
