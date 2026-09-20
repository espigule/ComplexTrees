#!/usr/bin/env python3
"""Finite degree certificates for nonembedded binary zippers.

Search is floating point and is NEVER itself a certificate. Verification uses
Gaussian integer polygon vertices and Fraction comparisons only. A successful
JSON certificate proves noninjectivity throughout its indicated parameter disk.
Words use symbols 0 and 1.  All dyadic source cylinders must be strictly disjoint.
"""
from __future__ import annotations
import argparse, json, math, cmath, collections
from fractions import Fraction as Q
from pathlib import Path
import numpy as np

def spec(name):
    if name not in {k+e for k in ('DD','DO','OO') for e in ('00','01','10','11')}:
        raise ValueError('Family must be one of the twelve marked DD/DO/OO presentations.')
    return {'DD':(0,0),'DO':(0,1),'OO':(1,1)}[name[:2]], tuple(map(int,name[2:]))

def maps(name,q):
    k,e=spec(name)
    return [(e[0]*q,(-1)**e[0]*q,k[0]),(q+e[1]*(1-q),(-1)**e[1]*(1-q),k[1])]

def compose(f,g):
    t,a,k=f; s,b,l=g
    return t+a*(s.conjugate() if k else s),a*(b.conjugate() if k else b),k^l

def wordmap(f,w):
    a=(0j,1+0j,0)
    for c in w:a=compose(a,f[int(c)])
    return a

def apply(f,z):
    t,a,k=f
    return t+a*(np.conjugate(z) if k else z)

def source_interval(name,w):
    _,e=spec(name); t=Q(0); a=Q(1)
    for c in w:
        i=int(c);t=t+a*Q(i+e[i],2);a=a*Q((-1)**e[i],2)
    return min(t,t+a),max(t,t+a)

def separated_intervals(name,u,v):
    I=source_interval(name,u);J=source_interval(name,v)
    return I[1]<J[0] or J[1]<I[0]

def curve(name,q,n):
    f=maps(name,q); _,e=spec(name);p=np.array([0j,1+0j])
    for _ in range(n):
        a=apply(f[0],p[::-1] if e[0] else p)
        b=apply(f[1],p[::-1] if e[1] else p)
        p=np.concatenate([a,b[1:]])
    return p

def boundary(a,b):
    return np.concatenate([a-b[0],(a[-1]-b)[1:],(a[::-1]-b[-1])[1:],(a[0]-b[::-1])[1:]])

def numerical_degree(p):
    a=p[:-1]; b=p[1:];d=b-a
    t=np.clip(-(a.real*d.real+a.imag*d.imag)/(np.abs(d)**2+1e-300),0,1)
    clearance=np.min(np.abs(a+t*d))
    angle=np.angle(b*np.conjugate(a)).sum()/(2*np.pi)
    return int(round(angle)),float(clearance)

def chord_cross(a,b,c,d):
    cross=lambda x,y:x.real*y.imag-x.imag*y.real
    return cross(b-a,c-a)*cross(b-a,d-a)<0 and cross(d-c,a-c)*cross(d-c,b-c)<0

def search(name,q,budget=30000,maxword=24,depths=(6,9,12)):
    f=maps(name,q);r=max(abs(q),abs(1-q))
    if r>=1:return None
    R=max(abs(1-q)/(2*(1-abs(q))),abs(q)/(2*(1-abs(1-q))))
    polys={n:curve(name,q,n) for n in depths}
    todo=collections.deque([('0','1')]);visits=0
    while todo and visits<budget:
        u,v=todo.popleft();visits+=1
        U=wordmap(f,u);V=wordmap(f,v)
        ru,rv=abs(U[1]),abs(V[1])
        if abs(apply(U,.5)-apply(V,.5))>R*(ru+rv):continue
        if separated_intervals(name,u,v):
            endsu=apply(U,np.array([0j,1+0j]));endsv=apply(V,np.array([0j,1+0j]))
            if chord_cross(*endsu,*endsv):
                for n,p in polys.items():
                    P=boundary(apply(U,p),apply(V,p));d,margin=numerical_degree(P)
                    tail=max(ru,rv)*r**n*abs(q-.5)/(1-r)
                    if d and margin>tail*1.1:
                        return dict(family=name,q=[q.real,q.imag],u=u,v=v,polygon_depth=n,numerical_degree=d,numerical_clearance=margin,numerical_tail=tail,visits=visits)
        if max(len(u),len(v))>=maxword:continue
        if ru>=rv:todo.extend([(u+'0',v),(u+'1',v)])
        else:todo.extend([(u,v+'0'),(u,v+'1')])
    return dict(family=name,q=[q.real,q.imag],status='unresolved',visits=visits,pending=len(todo))

# Exact Gaussian-integer arithmetic. Every polygon has a common denominator.
def gmul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gadd(a,b):return (a[0]+b[0],a[1]+b[1])
def gsub(a,b):return (a[0]-b[0],a[1]-b[1])
def gconj(a):return a[0],-a[1]

def exact_base(name,x,y):
    D=math.lcm(x.denominator,y.denominator);X=x.numerator*(D//x.denominator);Y=y.numerator*(D//y.denominator)
    k,e=spec(name); q=(X,Y);b=(D-X,-Y)
    f=[((e[0]*X,e[0]*Y),((-1)**e[0]*X,(-1)**e[0]*Y),k[0]),
       (gadd(q,(e[1]*b[0],e[1]*b[1])),((-1)**e[1]*b[0],(-1)**e[1]*b[1]),k[1])]
    return D,f

def exact_apply(f,p,den,D):
    t,a,k=f
    return [(gadd((t[0]*den,t[1]*den),gmul(a,gconj(z) if k else z))) for z in p],den*D

def exact_curve(name,x,y,n):
    D,f=exact_base(name,x,y);_,e=spec(name);p=[(0,0),(1,0)];den=1
    for _ in range(n):
        a,_=exact_apply(f[0],p[::-1] if e[0] else p,den,D)
        b,newden=exact_apply(f[1],p[::-1] if e[1] else p,den,D)
        assert a[-1]==b[0]
        p=a+b[1:];den=newden
    return p,den,D,f

def exact_word_polygon(f,w,p,den,D):
    for c in w[::-1]:p,den=exact_apply(f[int(c)],p,den,D)
    return p,den

def rational_upper_sqrt(q,den=10**12):
    # Exact ceil(sqrt(q)*den)/den.
    v=(q.numerator*den*den)//q.denominator;n=math.isqrt(v)
    while Q(n*n,den*den)<q:n+=1
    return Q(n,den)

def distance_squared_segment(a,b):
    d=gsub(b,a);dd=d[0]*d[0]+d[1]*d[1];dot=a[0]*d[0]+a[1]*d[1]
    if dd==0 or dot>=0:return Q(a[0]*a[0]+a[1]*a[1])
    if -dot>=dd:return Q(b[0]*b[0]+b[1]*b[1])
    cross=a[0]*d[1]-a[1]*d[0]
    return Q(cross*cross,dd)

def exact_winding(p):
    w=0
    for a,b in zip(p,p[1:]):
        cross=a[0]*b[1]-a[1]*b[0]
        if a[1]<=0<b[1] and cross>0:w+=1
        elif b[1]<=0<a[1] and cross<0:w-=1
    return w

def fraction_json(x):return str(x)

def verify(name,x,y,u,v,n,radius=Q(0)):
    spec(name)
    if radius<0 or n<0:raise ValueError('Radius and polygon depth must be nonnegative.')
    if not u or not v or set(u+v)-{'0','1'}:raise ValueError('Words must be nonempty binary strings.')
    assert separated_intervals(name,u,v), 'Source parameter cylinders must be strictly disjoint.'
    p,den,D,f=exact_curve(name,x,y,n)
    a,da=exact_word_polygon(f,u,p,den,D);b,db=exact_word_polygon(f,v,p,den,D)
    common=math.lcm(da,db);a=[(z[0]*(common//da),z[1]*(common//da)) for z in a];b=[(z[0]*(common//db),z[1]*(common//db)) for z in b]
    P=[gsub(z,b[0]) for z in a]+[gsub(a[-1],z) for z in b[1:]]+[gsub(z,b[-1]) for z in a[-2::-1]]+[gsub(a[0],z) for z in b[-2::-1]]
    assert P[0]==P[-1]
    r1=rational_upper_sqrt(x*x+y*y);r2=rational_upper_sqrt((1-x)**2+y*y);r=max(r1,r2)
    assert r+radius<1
    delta=rational_upper_sqrt((x-Q(1,2))**2+y*y)
    scale=lambda w:r1**w.count('0')*r2**w.count('1')
    tail=max(scale(u),scale(v))*r**n*delta/(1-r)
    R=max(r2/(2*(1-r1)),r1/(2*(1-r2)))
    lipschitz=(R+Q(1,2))/(1-r-radius)
    error=tail+2*lipschitz*radius
    error_scaled_sq=(error*common)**2
    minimum=None;passed=True
    for a,b in zip(P,P[1:]):
        ds=distance_squared_segment(a,b)
        if minimum is None or ds<minimum:minimum=ds
        if ds<=error_scaled_sq:passed=False
    degree=exact_winding(P)
    passed=passed and degree!=0
    I=source_interval(name,u);J=source_interval(name,v)
    return dict(schema='zipper-degree-capture-v1',status='proved_nonembedded_parameter_disk' if passed else 'unresolved',
        family=name,q={'real':str(x),'imag':str(y)},parameter_disk_radius=str(radius),
        words={'u':u,'v':v},source_intervals=[[str(t) for t in I],[str(t) for t in J]],
        polygon_depth=n,degree=degree,polygon_segments=len(P)-1,
        contraction_upper_bounds=[str(r1),str(r2)],curve_tail_bound=str(tail),parameter_motion_bound=str(2*lipschitz*radius),
        minimum_polygon_boundary_distance_squared=str(minimum/(common*common)),
        combined_error_squared=str(error*error),
        approximate_clearance=math.sqrt(float(minimum/(common*common))),approximate_combined_error=float(error),
        arithmetic='Exact Python integers and fractions.Fraction; square-root bounds rounded upwards by integer square comparisons.',
        conclusion='The canonical zipper map is not injective at every parameter in the indicated closed complex disk.' if passed else 'No conclusion.')

def main():
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='cmd',required=True)
    a=s.add_parser('search');a.add_argument('family');a.add_argument('real');a.add_argument('imag');a.add_argument('--budget',type=int,default=30000);a.add_argument('--maxword',type=int,default=24)
    b=s.add_parser('verify');b.add_argument('family');b.add_argument('real');b.add_argument('imag');b.add_argument('u');b.add_argument('v');b.add_argument('depth',type=int);b.add_argument('--radius',default='0')
    args=p.parse_args()
    if args.cmd=='search':o=search(args.family,complex(float(Q(args.real)),float(Q(args.imag))),args.budget,args.maxword)
    else:o=verify(args.family,Q(args.real),Q(args.imag),args.u,args.v,args.depth,Q(args.radius))
    print(json.dumps(o,indent=2))
if __name__=='__main__':main()
