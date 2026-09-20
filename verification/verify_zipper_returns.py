"""Exact integer-polynomial checks for every displayed zipper return identity.

Only Python standard library is required. Variables q and b are independent;
conjugation exchanges them, then b=conjugate(q) gives the planar identities.
A map is represented by c+a*J_k(z). These are entire parameter-family checks,
not numerical samples.
"""
from pathlib import Path
from fractions import Fraction
import json

class P:
    def __init__(self, terms):
        if isinstance(terms, P): terms=terms.t
        if isinstance(terms, (int, Fraction)): terms={(0,0):Fraction(terms)}
        self.t={k:Fraction(v) for k,v in terms.items() if v}
    def __add__(self, other):
        out=dict(self.t)
        for k,v in P(other).t.items():out[k]=out.get(k,0)+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.t.items()})
    def __sub__(self, other):return self+-P(other)
    def __rsub__(self, other):return P(other)+-self
    def __mul__(self, other):
        out={}
        for (i,j),a in self.t.items():
            for (k,l),b in P(other).t.items():
                out[(i+k,j+l)]=out.get((i+k,j+l),0)+a*b
        return P(out)
    __rmul__=__mul__
    def __truediv__(self, n):return self*Fraction(1,n)
    def conj(self):return P({(j,i):v for (i,j),v in self.t.items()})
    def __eq__(self, other):return self.t==P(other).t
    def __repr__(self):return str({str(k):str(v) for k,v in sorted(self.t.items())})

q=P({(1,0):1}); qb=q.conj()

def J(z,k):return z if k==1 else z.conj()

def compose(f,g):
    c,a,k=f;d,b,l=g
    return c+a*J(d,k),a*J(b,k),k*l

def maps(kind,signature):
    k={'DD':(1,1),'DO':(1,-1),'OO':(-1,-1)}[kind]
    e1,e2=map(int,signature)
    return ((e1*q,(-1)**e1*q,k[0]),
            (q+e2*(1-q),(-1)**e2*(1-q),k[1]))

def word(fs,s):
    f=(P(0),P(1),1)
    for c in s:f=compose(f,fs[int(c)-1])
    return f

rows=[
 ('DD','01','12','22','1','1',q),
 ('DD','10','11','21','2','2',1-q),
 ('DD','11','1','2','12','21',q*(1-q)),
 ('DO','10','11','21','22','22',(1-q)*(1-qb)),
 ('DO','11','1','2','1212','2121',q*qb*(1-q)*(1-qb)),
 ('OO','01','12','22','11','11',q*qb),
 ('OO','10','11','21','22','22',(1-q)*(1-qb))]

results=[]
for kind,e,u,v,r,s,eta in rows:
    fs=maps(kind,e);H=(q-eta*q,eta,1)
    for seed,tail in ((u,r),(v,s)):
        lhs=word(fs,seed+tail);rhs=compose(H,word(fs,seed))
        assert lhs==rhs, (kind,e,seed,tail,lhs,rhs)
    results.append(dict(plane=kind,signature=e,U=u,V=v,r=r,s=s,
                        identity='symbolically verified',eta_polynomial=repr(eta)))

for kind in ('DD','DO','OO'):
    for e in ('00','01','10','11'):
        f,g=maps(kind,e)
        assert f[0]+f[1]/2==q/2
        assert g[0]+g[1]/2==(1+q)/2
f,g=maps('OO','00')
assert f[0]+f[1]*qb==q*qb
assert g[0]+g[1]*qb==q+qb-q*qb
out=dict(status='passed',exact_return_identities=14,exact_centre_identities=24,
         exact_triangle_vertex_identities=2,rows=results)
Path(__file__).with_name('zipper_return_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
